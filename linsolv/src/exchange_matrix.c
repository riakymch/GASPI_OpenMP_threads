/*
 * exchange_matrix.c
 */

#include "exchange_matrix.h"

#include <mpi.h>
#include <GASPI.h>
#include <stddef.h>
#include "util.h"
#include "setup_comm.h"

static gaspi_offset_t *send_byteoffset;
static gaspi_offset_t *recv_byteoffset;
static gaspi_offset_t *remote_byteoffset;
static gaspi_notification_id_t *notification;
static gaspi_number_t nqueues;
static gaspi_number_t queue_size_max;
static gaspi_queue_id_t queue_id = 0;

#ifdef USE_GASPI_TEST
static int counter = 0;
#endif

#define NSEGM 2

/*******************************************************************************
*
*******************************************************************************/
void matrix_comm_init(CommMap *comap, int ncols)
{
  CHECK(comap != NULL);

  const int ndom = comap->ndomains;
  const int ncommdom = comap->ncommdomains;
  const int myid = comap->thisdomain;
  const int tag = 47;

  gaspi_number_t notif_num;

  GASPICHECK(gaspi_queue_num(&nqueues));
  GASPICHECK(gaspi_queue_size_max(&queue_size_max));
  CHECK(queue_size_max >= (gaspi_number_t) (2 * ncommdom));
  GASPICHECK(gaspi_notification_num(&notif_num));
  CHECK(notif_num >= (gaspi_number_t) ndom);

  send_byteoffset = check_malloc(ndom * sizeof(gaspi_offset_t));
  recv_byteoffset = check_malloc(ndom * sizeof(gaspi_offset_t));
  remote_byteoffset = check_malloc(ndom * sizeof(gaspi_offset_t));
  notification = check_malloc(ndom * sizeof(gaspi_offset_t));

  /* calculate local send and recv offsets */
  gaspi_offset_t nsend = 0;
  gaspi_offset_t nrecv = 0;
  for(int i = 0; i < ncommdom; i++)
  {
    const int k = comap->commpartner[i];

    send_byteoffset[k] = nsend;
    recv_byteoffset[k] = nrecv;
    nsend += comap->sendcount[k] * ncols * sizeof(double);
    nrecv += comap->recvcount[k] * ncols * sizeof(double);
  }

  /* create gaspi segments */
  for(int id = 0; id < NSEGM; id++)
  {
    GASPICHECK(gaspi_segment_create(id, (gaspi_size_t) nsend,
               GASPI_GROUP_ALL, GASPI_BLOCK, GASPI_ALLOC_DEFAULT));

    GASPICHECK(gaspi_segment_create(id + NSEGM, (gaspi_size_t) nrecv,
               GASPI_GROUP_ALL, GASPI_BLOCK, GASPI_ALLOC_DEFAULT));
  }

  /* exchange remote offsets */
  for(int i = 0; i < ncommdom; i++)
  {
    const int k = comap->commpartner[i];

    if(k > myid)
    {
      MPI_Send(&recv_byteoffset[k], sizeof(gaspi_offset_t), MPI_BYTE, k,
               tag, MPI_COMM_WORLD);
      MPI_Recv(&remote_byteoffset[k], sizeof(gaspi_offset_t), MPI_BYTE, k,
               tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
      MPI_Recv(&remote_byteoffset[k], sizeof(gaspi_offset_t), MPI_BYTE, k,
               tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      MPI_Send(&recv_byteoffset[k], sizeof(gaspi_offset_t), MPI_BYTE, k,
               tag, MPI_COMM_WORLD);
    }
  }

  /* exchange notification IDs */
  for(int i = 0; i < ncommdom; i++)
  {
    const int k = comap->commpartner[i];
    gaspi_notification_id_t notif_id = i;

    if(k > myid)
    {
      MPI_Send(&notif_id, sizeof(gaspi_notification_id_t), MPI_BYTE, k,
               tag, MPI_COMM_WORLD);
      MPI_Recv(&notification[k], sizeof(gaspi_notification_id_t), MPI_BYTE, k,
               tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
      MPI_Recv(&notification[k], sizeof(gaspi_notification_id_t), MPI_BYTE, k,
               tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      MPI_Send(&notif_id, sizeof(gaspi_notification_id_t), MPI_BYTE, k,
               tag, MPI_COMM_WORLD);
    }
  }

#ifdef USE_GASPI_TEST
  /* check number of recv domains matches number of comm domains */
  for(int i = 0; i < ncommdom; i++)
    CHECK(comap->recvcount[comap->commpartner[i]] > 0);
#endif

} /* matrix_comm_init() */

/*******************************************************************************
*
*******************************************************************************/
void matrix_comm_end(void)
{
  check_free(send_byteoffset);
  check_free(recv_byteoffset);
  check_free(remote_byteoffset);
  check_free(notification);

  for(int id = 0; id < NSEGM; id++)
  {
    GASPICHECK(gaspi_segment_delete(id));
    GASPICHECK(gaspi_segment_delete(id + NSEGM));
  }
} /* matrix_comm_end() */

/*******************************************************************************
*
*******************************************************************************/
void exchange_matrix(const CommMap *comap, Matrix matrix, int stage)
{
  if(comap == NULL || comap->ndomains <= 1)
    return;

  const int ncols = matrix.cols;
  double **mat = matrix.m;
  const int ncommpartner = comap->ncommdomains;
  int *commpartner = comap->commpartner;

  const int send_segm_id = stage % NSEGM;
  const int recv_segm_id = send_segm_id + NSEGM;

#ifdef USE_GASPI_TEST
# pragma omp single nowait
  counter = 0;
#endif

# pragma omp single
  {
    gaspi_number_t queue_size;
    GASPICHECK(gaspi_queue_size(queue_id, &queue_size));
    if(queue_size + 2 * ncommpartner > queue_size_max)
    {
      queue_id = (queue_id + 1) % nqueues;
      GASPICHECK(gaspi_wait(queue_id, GASPI_BLOCK));
    }
  }

  /* copy in and send data */
# pragma omp for COMM_SCHEDULE nowait
  for(int i = 0; i < ncommpartner; i++)
  {
    const int dest = commpartner[i];
    const int sendcount = comap->sendcount[dest];

    if(sendcount > 0)
    {
      gaspi_pointer_t segm_ptr;
      GASPICHECK(gaspi_segment_ptr(send_segm_id, &segm_ptr));
      double *sbuf = (double *) ((char *) segm_ptr + send_byteoffset[dest]);
      int *sendindex = comap->sendindex[dest];

      /* copy data to send segment */
      for(int j = 0; j < sendcount; j++)
        for(int col = 0; col < ncols; col++)
          sbuf[j * ncols + col] = mat[sendindex[j]][col];

      gaspi_size_t size = sendcount * ncols * sizeof(double);

      GASPICHECK(gaspi_write_notify(send_segm_id, send_byteoffset[dest], dest,
                                    recv_segm_id, remote_byteoffset[dest],
                                    size, notification[dest], 1, queue_id,
                                    GASPI_BLOCK));
    }
  }

  /* wait for data arrival and copy out */
#ifndef USE_GASPI_TEST /* wait on single notification using GASPI_BLOCK */
# pragma omp for COMM_SCHEDULE
  for(int i = 0; i < ncommpartner; i++)
  {
    const int source = commpartner[i];
    const int recvcount = comap->recvcount[source];

    if(recvcount > 0)
    {
      gaspi_notification_id_t notif_id;
      gaspi_notification_t notif_val;

      GASPICHECK(gaspi_notify_waitsome(recv_segm_id, i, 1, &notif_id,
                                       GASPI_BLOCK));
      GASPICHECK(gaspi_notify_reset(recv_segm_id, notif_id, &notif_val));

      gaspi_pointer_t segm_ptr;
      GASPICHECK(gaspi_segment_ptr(recv_segm_id, &segm_ptr));
      double *rbuf = (double *) ((char *) segm_ptr + recv_byteoffset[source]);
      int *recvindex = comap->recvindex[source];

      /* copy data from recv segment */
      for(int j = 0; j < recvcount; j++)
        for(int col = 0; col < ncols; col++)
          mat[recvindex[j]][col] = rbuf[j * ncols + col];
    }
  }
#else /* wait on notification range using GASPI_TEST */
  int local_counter;

# pragma omp atomic read
  local_counter = counter;

  while(local_counter < ncommpartner)
  {
    gaspi_notification_id_t notif_id;
    gaspi_notification_t notif_val = 0;

    gaspi_return_t ret = gaspi_notify_waitsome(recv_segm_id, 0, ncommpartner,
                                               &notif_id, GASPI_TEST);

    if(ret != GASPI_TIMEOUT)
    {
      CHECK(ret == GASPI_SUCCESS);
      GASPICHECK(gaspi_notify_reset(recv_segm_id, notif_id, &notif_val));
    }

    if(notif_val)
    {
#     pragma omp atomic
      counter++;

      const int source = commpartner[notif_id];
      const int recvcount = comap->recvcount[source];

      if(recvcount > 0)
      {
        gaspi_pointer_t segm_ptr;
        GASPICHECK(gaspi_segment_ptr(recv_segm_id, &segm_ptr));
        double *rbuf = (double *) ((char *) segm_ptr + recv_byteoffset[source]);
        int *recvindex = comap->recvindex[source];

        /* copy data from recv segment */
        for(int j = 0; j < recvcount; j++)
          for(int col = 0; col < ncols; col++)
            mat[recvindex[j]][col] = rbuf[j * ncols + col];
      }
    }
#   pragma omp atomic read
    local_counter = counter;
  }
# pragma omp barrier
#endif
} /* exchange_matrix() */
