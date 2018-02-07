#include <iostream>
#include <iomanip>
#include <cassert>
#include <thread>
#include <stdexcept>
#include <algorithm>
#include <string>
#include "mpi_interface.h"

using namespace std;

MPI_Interface::MPI_Interface()
 : cnt_(0), mpimutex_() {

#ifdef HAVE_MPI_H
  int provided;
  MPI_Init_thread(0, 0, MPI_THREAD_MULTIPLE, &provided);
  if (provided != MPI_THREAD_MULTIPLE)
    throw runtime_error("MPI_THREAD_MULTIPLE not provided");

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size_);
  rank_ = world_rank_;
  size_ = world_size_;

  // print out the node name
  {
    constexpr const size_t maxlen = MPI_MAX_PROCESSOR_NAME;
    int len;
    char name[maxlen];
    MPI_Get_processor_name(name, &len);

    unique_ptr<char[]> buf(new char[maxlen*size_]);
    unique_ptr<int[]> lens(new int[size_]);
    MPI_Gather(static_cast<void*>(name), maxlen, MPI_CHAR, static_cast<void*>(buf.get()), maxlen, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Gather(static_cast<void*>(&len),      1, MPI_INT,  static_cast<void*>(lens.get()),     1, MPI_INT,  0, MPI_COMM_WORLD);
    if (rank() == 0) {
      for (int i = 0; i != size_; ++i)
        cout << left << "    " << setw(32) << string(&buf[i*maxlen], &buf[i*maxlen+lens[i]]) << right << endl;
      cout << endl;
    }
  }

  // obtain the upper bound of tags
  {
    int flag, *get_val;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &get_val, &flag);
    assert(flag && *get_val >= 32767); // this is what the standard says
    tag_ub_ = *get_val;
  }

  // set MPI_COMM_WORLD to mpi_comm_
  mpi_comm_ = MPI_COMM_WORLD;
  
  //Printing
  if (rank() != 0) {
    cout_orig = cout.rdbuf();
    cout.rdbuf(ss_.rdbuf());
  }
#else
  rank_ = 0;
  size_ = 1;
#endif
}


MPI_Interface::~MPI_Interface() {
#ifdef HAVE_MPI_H
  MPI_Finalize();
#endif
}


void MPI_Interface::barrier() const {
#ifdef HAVE_MPI_H
  MPI_Barrier(mpi_comm_);
#endif
}


void MPI_Interface::allreduce(double* a, const size_t size) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_DOUBLE, MPI_SUM, mpi_comm_);
#endif
}


void MPI_Interface::allreduce(int* a, const size_t size) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_INT, MPI_SUM, mpi_comm_);
#endif
}


void MPI_Interface::allreduce(complex<double>* a, const size_t size) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_);
#endif
}

void MPI_Interface::broadcast(int* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  assert(size != 0); 
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Bcast(static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_INT, root, mpi_comm_);
#endif
}

void MPI_Interface::broadcast(size_t* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  static_assert(sizeof(size_t) == sizeof(unsigned long long), "size_t is assumed to be the same size as unsigned long long");
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Bcast(static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_UNSIGNED_LONG_LONG, root, mpi_comm_);
#endif
}

void MPI_Interface::broadcast(double* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Bcast(static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_DOUBLE, root, mpi_comm_);
#endif
}


void MPI_Interface::broadcast(complex<double>* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Bcast(static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_CXX_DOUBLE_COMPLEX, root, mpi_comm_);
#endif
}


void MPI_Interface::allgather(const double* send, const size_t ssize, double* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_DOUBLE, static_cast<void*>(rec), rsize, MPI_DOUBLE, mpi_comm_);
#else
  assert(ssize == rsize);
  copy_n(send, ssize, rec);
#endif
}


void MPI_Interface::allgather(const complex<double>* send, const size_t ssize, complex<double>* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_CXX_DOUBLE_COMPLEX, static_cast<void*>(rec), rsize, MPI_CXX_DOUBLE_COMPLEX, mpi_comm_);
#else
  assert(ssize == rsize);
  copy_n(send, ssize, rec);
#endif
}


void MPI_Interface::allgather(const size_t* send, const size_t ssize, size_t* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  static_assert(sizeof(size_t) == sizeof(unsigned long long), "size_t is assumed to be the same size as unsigned long long");
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_UNSIGNED_LONG_LONG, static_cast<void*>(rec), rsize, MPI_UNSIGNED_LONG_LONG, mpi_comm_);
#else
  assert(ssize == rsize);
  copy_n(send, ssize, rec);
#endif
}


void MPI_Interface::allgather(const int* send, const size_t ssize, int* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_INT, static_cast<void*>(rec), rsize, MPI_INT, mpi_comm_);
#else
  assert(ssize == rsize);
  copy_n(send, ssize, rec);
#endif
}

void MPI_Interface::ibroadcast(double* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    MPI_Ibcast(static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_DOUBLE, root, mpi_comm_, &c);
  }
#endif
}

int MPI_Interface::request_send(const double* sbuf, const size_t size, const int dest, const int tag) {
#ifdef HAVE_MPI_H
  assert(tag <= tag_ub_);
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    // I hate const_cast. Blame the MPI C binding
    MPI_Isend(const_cast<double*>(sbuf+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_DOUBLE, dest, tag, mpi_comm_, &c);
    rq.push_back(c);
  }
#endif
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  request_.emplace(cnt_, rq);
#endif
  ++cnt_;
  return cnt_-1;
}


int MPI_Interface::request_send(const complex<double>* sbuf, const size_t size, const int dest, const int tag) {
#ifdef HAVE_MPI_H
  assert(tag <= tag_ub_);
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    // I hate const_cast. Blame the MPI C binding
    MPI_Isend(const_cast<complex<double>*>(sbuf+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_CXX_DOUBLE_COMPLEX, dest, tag, mpi_comm_, &c);
    rq.push_back(c);
  }
#endif
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  request_.emplace(cnt_, rq);
#endif
  ++cnt_;
  return cnt_-1;
}


int MPI_Interface::request_send(const size_t* sbuf, const size_t size, const int dest, const int tag) {
#ifdef HAVE_MPI_H
  assert(tag <= tag_ub_);
  static_assert(sizeof(size_t) == sizeof(unsigned long long), "size_t is assumed to be the same size as unsigned long long");
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    // I hate const_cast. Blame the MPI C binding
    MPI_Isend(const_cast<size_t*>(sbuf+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_UNSIGNED_LONG_LONG, dest, tag, mpi_comm_, &c);
    rq.push_back(c);
  }
#endif
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  request_.emplace(cnt_, rq);
#endif
  ++cnt_;
  return cnt_-1;
}



int MPI_Interface::request_recv(double* rbuf, const size_t size, const int origin, const int tag) {
#ifdef HAVE_MPI_H
  assert(tag <= tag_ub_);
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    MPI_Irecv(rbuf+i*bsize, (i+1 == nbatch ? size-i*bsize : bsize), MPI_DOUBLE, (origin == -1 ? MPI_ANY_SOURCE : origin), (tag==-1 ? MPI_ANY_TAG : tag), mpi_comm_, &c);
    rq.push_back(c);
  }
#endif
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  request_.emplace(cnt_, rq);
#endif
  ++cnt_;
  return cnt_-1;
}


int MPI_Interface::request_recv(complex<double>* rbuf, const size_t size, const int origin, const int tag) {
#ifdef HAVE_MPI_H
  assert(tag <= tag_ub_);
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    MPI_Irecv(rbuf+i*bsize, (i+1 == nbatch ? size-i*bsize : bsize), MPI_CXX_DOUBLE_COMPLEX, (origin == -1 ? MPI_ANY_SOURCE : origin), (tag==-1 ? MPI_ANY_TAG : tag), mpi_comm_, &c);
    rq.push_back(c);
  }
#endif
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  request_.emplace(cnt_, rq);
#endif
  ++cnt_;
  return cnt_-1;
}


int MPI_Interface::request_recv(size_t* rbuf, const size_t size, const int origin, const int tag) {
#ifdef HAVE_MPI_H
  assert(tag <= tag_ub_);
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    MPI_Irecv(rbuf+i*bsize, (i+1 == nbatch ? size-i*bsize : bsize), MPI_UNSIGNED_LONG_LONG, (origin == -1 ? MPI_ANY_SOURCE : origin), (tag==-1 ? MPI_ANY_TAG : tag), mpi_comm_, &c);
    rq.push_back(c);
  }
#endif
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  request_.emplace(cnt_, rq);
#endif
  ++cnt_;
  return cnt_-1;
}


void MPI_Interface::wait(const int rq) {
#ifdef HAVE_MPI_H
  lock_guard<mutex> lock(mpimutex_);
  auto i = request_.find(rq);
  assert(i != request_.end());
  for (auto& j : i->second)
    MPI_Wait(&j, MPI_STATUS_IGNORE);
#endif
}


void MPI_Interface::cancel(const int rq) {
#ifdef HAVE_MPI_H
  lock_guard<mutex> lock(mpimutex_);
  auto i = request_.find(rq);
  assert(i != request_.end());
  for (auto& j : i->second)
    MPI_Cancel(&j);
#endif
}


bool MPI_Interface::test(const int rq) {
  bool out = true;
#ifdef HAVE_MPI_H
  lock_guard<mutex> lock(mpimutex_);
  auto i = request_.find(rq);
  assert(i != request_.end());
  for (auto& j : i->second) {
    int b;
    MPI_Test(&j, &b, MPI_STATUS_IGNORE);
    out &= b;
  }
#endif
  return out;
}


// MPI Communicators
void MPI_Interface::split(const int n) {
#ifdef HAVE_MPI_H
  MPI_Comm new_comm;
  const int icomm = world_rank_ / n;
  MPI_Comm_split(MPI_COMM_WORLD, icomm, world_rank_, &new_comm);
  mpi_comm_ = new_comm;
  MPI_Comm_rank(mpi_comm_, &rank_);
  MPI_Comm_size(mpi_comm_, &size_);
#endif
}


void MPI_Interface::merge() {
#ifdef HAVE_MPI_H
  MPI_Comm_free(&mpi_comm_);
  mpi_comm_ = MPI_COMM_WORLD;
  rank_ = world_rank_;
  size_ = world_size_;
#endif
}
