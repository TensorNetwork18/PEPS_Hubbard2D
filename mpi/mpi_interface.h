#ifndef _MPI_INTERFACE_H_
#define _MPI_INTERFACE_H_

#include <stddef.h>
#include <memory>
#include <complex>
#include <mutex>
#include <vector>
#include <map>
#ifdef HAVE_MPI_H
 #include <mpi.h>
#endif


class MPI_Interface {
  protected:
    int world_rank_;
    int world_size_;
    int rank_;
    int size_;

    int cnt_;

    //priting
    std::streambuf* cout_orig;
    std::stringstream ss_;
    // request handles
#ifdef HAVE_MPI_H
    MPI_Comm mpi_comm_;
    std::map<int, std::vector<MPI_Request>> request_;
#endif

    // maximum size of the MPI buffer
    static constexpr size_t bsize = 100000000LU;

    // mutex for isend and irecv
    mutable std::mutex mpimutex_;

    // MPI's internal variables
    int tag_ub_;

  public:
    MPI_Interface();
    ~MPI_Interface();

    int world_rank() const { return world_rank_; }
    int world_size() const { return world_size_; }
    int rank() const { return rank_; }
    int size() const { return size_; }
    bool last() const { return rank() == size()-1; }

    // collective functions
    // barrier
    void barrier() const;
    // sum reduce and broadcast to each process
    void allreduce(int*, const size_t size) const;
    void allreduce(double*, const size_t size) const;
    void allreduce(std::complex<double>*, const size_t size) const;
    // broadcast
    void broadcast(int*, const size_t size, const int root) const;
    void broadcast(size_t*, const size_t size, const int root) const;
    void broadcast(double*, const size_t size, const int root) const;
    void broadcast(std::complex<double>*, const size_t size, const int root) const;
    void allgather(const double* send, const size_t ssize, double* rec, const size_t rsize) const;
    void allgather(const std::complex<double>* send, const size_t ssize, std::complex<double>* rec, const size_t rsize) const;
    void allgather(const size_t* send, const size_t ssize, size_t* rec, const size_t rsize) const;
    void allgather(const int* send, const size_t ssize, int* rec, const size_t rsize) const;

    void ibroadcast(double*, const size_t size, const int root) const;
    // one-sided communication with Isend, Irecv
    int request_send(const double* sbuf, const size_t size, const int dest, const int tag);
    int request_send(const std::complex<double>* sbuf, const size_t size, const int dest, const int tag);
    int request_send(const size_t* sbuf, const size_t size, const int dest, const int tag);
    int request_recv(double* rbuf, const size_t size, const int source = -1, const int tag = -1);
    int request_recv(std::complex<double>* rbuf, const size_t size, const int source = -1, const int tag = -1);
    int request_recv(size_t* rbuf, const size_t size, const int source = -1, const int tag = -1);
    void wait(const int rq);
    void cancel(const int rq);
    bool test(const int rq);

#ifdef HAVE_MPI_H
    // communicators. n is the number of processes per communicator.
    const MPI_Comm& mpi_comm() const { return mpi_comm_; } 
#endif
    void split(const int n);
    void merge();
};

extern MPI_Interface* mpi__;
#endif
