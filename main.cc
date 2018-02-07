#include <thread>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cassert>
#include <vector>
#include "static.h"
#include "sparse_tensor/qnum.h"
#include "sparse_tensor/basis.h"
#include "peps/peps.h"
#include "peps/gauge.h"
#include "hubbard/minimizeE.h"
#ifdef HAVE_TBB
  #include "thread/tbb_interface.h"
#endif

using namespace std;

int main(int argc, char** argv) {
  static_variables();
#ifdef HAVE_TBB
  tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic);
  //tbb::task_scheduler_init init(1);
#endif

  // INPUT //
  string inputfile;
  for (int i = 1; i != argc; ++i) {
    if (i+1 != argc) {
      if (string(argv[i]) == "-i") {
        inputfile = argv[i+1];
      }
    }
  }
  if (inputfile.empty())
    throw runtime_error("input file was not specified");

  // Parameters //
  bool do_minimize = false;
  string rfile, ifile;
  size_t Lx, Ly;
  double u_term;
  vector<Basis> basisX, basisY;
  map<Qnum, size_t> cMap;

  double gauge_thresh;

  double perturb;
  double thresh;
  size_t macro_cycle;
  size_t max_cycle;

  //READ INPUT
  ifstream input(inputfile);
  if (!input.good())
    throw runtime_error("could not open the input file " + inputfile);

  cout << "READ INPUT: " << inputfile << endl;
  string line;
  while (getline(input, line)) {
    istringstream iss(line);
    string key;
    if (getline(iss, key)) {

      if (key[0] == 'L') {

        cout << "Lattice" << endl;
        getline(input, line);
        iss.clear();
        iss.str(line);
        iss >> Lx >> Ly >> u_term;
        cout << Lx << " " << Ly << " " << u_term << endl << endl;

      } else if (key[0] == 'X') {

        cout << "XBasis" << endl;
        for (size_t i = 0; i != Lx*(Ly+1); ++i) {
          getline(input, line);
          iss.clear();
          iss.str(line);

          int qn;
          vector<Qnum> qnum;
          while (iss >> qn) {
            cout << qn << " ";
            qnum.push_back(Qnum(qn, true));
          } cout << endl << " ";
          if (qnum.size() == 0)
            throw runtime_error("Input XBasis is not valid");

          basisX.push_back(Basis(Basis::OUT, qnum)); 
        }
        cout << endl;

      } else if (key[0] == 'Y') {

        cout << "YBasis" << endl;
        for (size_t i = 0; i != Ly*(Lx+1); ++i) {
          getline(input, line);
          iss.clear();
          iss.str(line);

          int qn; 
          vector<Qnum> qnum;
          while (iss >> qn) {
            cout << qn << " " ; 
            qnum.push_back(Qnum(qn, true));
          } cout << endl;

          if (qnum.size() == 0)
            throw runtime_error("Input YBasis is not valid");

          basisY.push_back(Basis(Basis::OUT, qnum));
        }
        cout << endl;

      } else if (key[0] == 'C') {

        cout << "Compress" << endl;
        getline(input, line);
        iss.clear();
        iss.str(line);
        string qstring;
        int qn;
        vector<Qnum> qnum;
        while (iss >> qn) {
          cout << qn << " ";
          qnum.push_back(Qnum(qn, true));
        } cout << endl;

        getline(input, line);
        iss.clear();
        iss.str(line);
        size_t qd;
        vector<size_t> qdim;
        while (iss >> qd) {
          cout << qd << " ";
          qdim.push_back(qd);
        } cout << endl;

        if (qnum.size() != qdim.size())
          throw runtime_error("Input Compress is not valid");

        for (size_t i = 0; i != qnum.size(); ++i) {
          cMap.emplace(qnum[i], qdim[i]);
        }
        cout << endl;

      } else if (key[0] == 'G') {

        cout << "Gauge" << endl;
        getline(input, line);
        iss.clear();
        iss.str(line);
        iss >> gauge_thresh;
        cout << gauge_thresh << endl << endl;

      } else if (key[0] == 'M') {
        
        do_minimize = true;
        cout << "Minimize" << endl;
        getline(input, line);
        iss.clear();
        iss.str(line);
        iss >> perturb >> thresh >> macro_cycle >> max_cycle;
        cout << perturb << " " << thresh << " " << macro_cycle << " " << max_cycle << endl << endl;

      } else if (key[0] == 'R') {

        getline(input, rfile);
        cout << "Restart" << endl;
        cout << rfile << endl << endl;

      } else if (key[0] == 'I') {
        
        getline(input, ifile);
        cout << "Init" << endl;
        cout << ifile << endl << endl;
      }
    }
  }

  shared_ptr<PEPS> ket;
  {
    //PEPS
    vector<Qtensor> Qts;
    Qnum q0(0, true);
    Qnum q1(1, true);
    Qnum q2(2, true);
    Basis s(Basis::OUT,{q0,q1,q1,q2});
    for (size_t ix = 0; ix != Lx; ++ix) {
      for (size_t iy = 0; iy != Ly; ++iy) {
        const size_t ax = ix * (Ly+1) + iy;
        const size_t ay = ix + iy * (Lx+1);
        Qts.push_back(Qtensor({ basisX[ax], basisY[ay], s, basisY[ay+1].reverse(), basisX[ax+1].reverse()}, true));
      }
    }

    ket = make_shared<PEPS>(Lx, Ly, Qts);
  }
 
  if (!ket)
   throw runtime_error("Input PEPS is not valid"); 

  if (!rfile.empty()) {

    cout << "PEPS restarted from : " << rfile << endl << endl; 
    ket->randomize(-0.05, 0.05);
    ket->load(rfile);

  } else if (!ifile.empty()) {
 
    cout << "PEPS initialized from : " << ifile << endl << endl;
    ket->randomize(-0.1, 0.1);   
    ifstream init(ifile);
    if (!init.good())
      throw runtime_error("could not open the INIT file " + ifile);
    while (getline(init, line)) {
      
      size_t sx, sy;
      int q0, q1, q2, q3, q4;
      size_t i0, i1, i2, i3, i4;
      double value;
      istringstream iss(line);
      iss.clear();
      iss.str(line);
      
      iss >> sx >> sy >> q0 >> q1 >> q2 >> q3 >> q4 >> i0 >> i1 >> i2 >> i3 >> i4 >> value;
      auto tensor_ptr = ket->site(sx, sy).get_block({Qnum(q0, true), Qnum(q1, true), Qnum(q2, true), Qnum(q3, true), Qnum(q4, true)});
      if (!tensor_ptr)
        throw runtime_error("Bad INIT file: " + line);
      tensor_ptr->element(i0,i1,i2,i3,i4) = value;
    }
     
  } else {
    
    cout << "PEPS randomized initialization " << endl << endl;
    ket->randomize(-1.0, 1.0);

  }

  if (do_minimize) {
    //OPTIMIZE//
    if (rfile.empty()) {  //I think it'd be better this way
      Hubbard2::MinimizeE obj0(u_term, ket, cMap);
      obj0.iteration(0.0, thresh, 1); 
      ket = obj0.output();
      tag__->reset();
    }

    for (size_t icycle = 0; icycle < max_cycle; icycle = icycle + macro_cycle) {
      if (Lx >1) {
        Gauge g(ket);
        g.iteration(gauge_thresh);
        ket = g.output();
      }    

      Hubbard2::MinimizeE obj(u_term, ket, cMap);
      auto chk = obj.iteration(perturb, thresh, macro_cycle);
      perturb = chk.second;
      ket = obj.output();
      tag__->reset();
      if (chk.first)
        break;
    }
  }

  //COMPUTE//
  if (Lx >1) {
    Gauge g(ket);
    g.iteration(gauge_thresh);
    ket = g.output();
  }

  Hubbard2::MinimizeE obj_last(u_term, ket, cMap);
  obj_last.compute();
  
  return 0;
}
