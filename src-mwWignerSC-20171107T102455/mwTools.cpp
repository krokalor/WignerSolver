// DO NOT edit outside SYNC/Sci/src/mwTools/ !

// Copyright 2008-2015 Maciej Woloszyn <woloszyn@agh.edu.pl>

// TODO tests for EACH function!!!

#include "mwTools.h"

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <stdexcept>
#include <unistd.h>
//#include <climits>
//#include <limits>

namespace mw {

////////////////////////////////////////////////////////////

std::string time()
{
  time_t curtime;
  struct tm *loctime;
  curtime = std::time (NULL);// Get the current time.
  loctime = std::localtime (&curtime); // Convert it to local time representation
  char date[MW_MAX_DATE_LEN]; // Buffer for date in format "%Y-%m-%d,%H:%M:%S"
  std::strftime(date,MW_MAX_DATE_LEN,"%Y-%m-%d,%H:%M:%S",loctime); // See strftime(3)
  std::string sTime(date);
  return sTime;
}
////////////////////////////////////////////////////////////

std::string host()
{
  char cstr_tmp[MW_CSTR_TMP_LEN];
  gethostname( cstr_tmp, MW_CSTR_TMP_LEN );
  std::string s( cstr_tmp );
  return s;
}

////////////////////////////////////////////////////////////

std::string commandLine(int argc, char *argv[])
{
  std::stringstream s;
  for(int i=0; i<argc; ++i)
    s << argv[i] << " ";
  return s.str();
}

////////////////////////////////////////////////////////////

std::string startTime()
{
  std::stringstream s;
  s << "## Start  (" << mw::host() << ") " << mw::time() ;
  return s.str();
}

////////////////////////////////////////////////////////////

std::string stopTime()
{
  std::stringstream s;
  s << "## Stop  (" << mw::host() << ") " << mw::time() ;
  return s.str();
}
////////////////////////////////////////////////////////////

std::string doubleToString(double x)
{
  std::stringstream s;
  s << x;
  return s.str();
}

////////////////////////////////////////////////////////////

v2D transpose(const v2D& data) {
  // this assumes that all inner vectors have the same size and
  // allocates space for the complete result in advance
  v2D result(data[0].size(),  v1D(data.size()));
  for (v1D::size_type i = 0; i < data[0].size(); i++) {
    for (v1D::size_type j = 0; j < data.size(); j++) {
      result[i][j] = data[j][i];
    }
  }
  return result;
}

    ////////////////////////////////////////////////////////////

std::string loadData(std::istream& in, std::vector<double>& v)
{
  v.clear();
  std::string info, line; //, sep;
  double x;
  getline(in, line);
  if(line.size()>0 && line[0]=='#') {
    info = line;
    getline(in, line);
  }
  do {
    std::istringstream iss(line);
    if(iss>>x) {
      v.push_back(x);
    }
  } while( getline(in, line) );
  return info;
}

////////////////////////////////////////////////////////////

std::string loadData(std::istream& in, std::vector<double>& vx,
                       std::vector<double>& vy)
{
  vx.clear();
  vy.clear();
  std::string info, line; //, sep;
  double x,y;
  getline(in, line);
  if(line.size()>0 && line[0]=='#') {
    info = line;
    getline(in, line);
  }
  do {
    std::istringstream iss(line);
    if(iss>>x>>y) {
      vx.push_back(x);
      vy.push_back(y);
    }
  } while( getline(in, line) );
  return info;
}

////////////////////////////////////////////////////////////

std::string loadData(std::istream& in, std::vector< std::vector<double> >& t)
{
  t.clear();
  std::string info, line; //, sep;
  double x;
  getline(in, line);
  if(line.size()>0 && line.at(0)=='#') {
    info = line;
    getline(in, line);
  }
  do {
    std::istringstream iss(line);
    if(iss >> x) {
      std::vector<double> row;
      row.push_back(x);
      while(iss >> x) {
        row.push_back(x);
      }
      t.push_back(row);
    }
  } while( getline(in, line) );
  return info;
}

////////////////////////////////////////////////////////////

void saveData(std::ostream &out, const std::vector<std::vector<double> >& v,
                std::string info)
{
    if(info.size() > 0 ) out << "# " << info << std::endl;

    const size_t Nx = v.size();
    for (size_t i=0; i!=Nx; ++i) {
        const size_t Ny = v.at(i).size();
        for (size_t j=0; j!=Ny; ++j) {
            out << v.at(i).at(j) << " ";
        }
        out << std::endl;
    }
}

////////////////////////////////////////////////////////////

void saveData(std::ostream& out, const std::vector<double>& vx,
                const std::vector<double>& vy, std::string info)
{
    size_t N = vx.size();
    if(N>vy.size()) N = vy.size();

    if(info.size() > 0 ) out << "# " << info << std::endl;

    for (size_t i=0; i!=N; ++i) {
        out << vx.at(i) << " " << vy.at(i) << std::endl;
    }
}

////////////////////////////////////////////////////////////

void saveData(std::ostream& out, const std::vector<double>& vx,
                std::string info)
{
    const size_t N = vx.size();

    if(info.size() > 0 ) out << "# " << info << std::endl;

    for (size_t i=0; i!=N; ++i) {
        out << vx.at(i) << std::endl;
    }
}

////////////////////////////////////////////////////////////

double linearSplineValue(double x, const std::vector<double>& xnodes,
                         const std::vector<double>& ynodes)
{
  size_t size = xnodes.size();
  double a,b;
// TODO move to new function linearSplineValid(...)
//   if(size > y.size()) {
//     size = y.size();
//     std::cerr << "WARNING: linearSpline: x.size() != y.size(), smaller ane will be taken.\n";
//   }
  //if(not sorted)...
  //double xp = std::numeric_limits<double>::min(); // for the previous x
  //if(x<xp)
  //  throw std::logic_error(std::string("load2dData: 1st column data not sorted!"));
  //xp=x;
  if(x <= xnodes[0])
    return ynodes[0];
  if(x >= xnodes[size-1])
    return ynodes[size-1];
  for(size_t i=0;i<(size-1);i++) {
    if(x<xnodes[i+1]) {
      a=(ynodes[i]-ynodes[i+1])/(xnodes[i]-xnodes[i+1]);
      b=(xnodes[i]*ynodes[i+1]-ynodes[i]*xnodes[i+1])/(xnodes[i]-xnodes[i+1]);
      return (a*x+b);
    }
  }
  // we should not never be here...
  throw std::logic_error(std::string("linearSplineValue: error"));
  return 0.0;
}

////////////////////////////////////////////////////////////

void discreteLinearSpline(const std::vector<double>& xnodes,
                          const std::vector<double>& ynodes, double dx,
                          std::vector<double>& vy)
{
  size_t size = xnodes.size();
  if(size < 1)
    throw std::logic_error(std::string(
      "discreteLinearSpline : invalid data on input"));
  double x0 = xnodes[0];
  double x1 = xnodes[size-1];
  vy.clear();
  size_t vy_estimated_size = size_t ((x1-x0)/dx) + 1;
  vy.reserve(vy_estimated_size);
  for(double x=x0; x<=x1; x+=dx) {
    vy.push_back( linearSplineValue(x, xnodes, ynodes) );
  }
}

////////////////////////////////////////////////////////////

std::vector<double> linearFit(const std::vector<double>& x,
        const std::vector<double>& y)
{
  const size_t N = x.size();
  if(N < 2)
    throw std::logic_error(std::string("linearFit: N < 2"));
  if(N != y.size())
    throw std::logic_error(std::string("linearFit: x.size() != y.size()"));
  const double xa = average(x);
  const double ya = average(y);
  double D = 0.;
  for(size_t i = 0; i < N; ++i) {
    D += (x[i]-xa)*(x[i]-xa);
  }
  double a = 0.;
  for(size_t i = 0; i < N; ++i) {
    a += y[i]*(x[i]-xa);
  }
  a /= D;
  double b = ya - a*xa;
  double Sy = 0.;
  for(size_t i = 0; i < N; ++i) {
    Sy += (y[i]-(a*x[i]+b))*(y[i]-(a*x[i]+b));
  }
  Sy = std::sqrt(Sy/(double)(N-2));
  double ua = Sy/std::sqrt(D);
  double ub = Sy*std::sqrt(1./(double)N + xa*xa/D);
  std::vector<double> result(4, 0);
  result[0] = a;
  result[1] = b;
  result[2] = ua;
  result[3] = ub;
  return result;
}

////////////////////////////////////////////////////////////

v1D derivative(const v1D& data, double dx, int) {
  v1D d(data.size()-1, 0.); // autocorr values
  for(size_t i = 0; i<d.size(); ++i) {
    d[i] = (data[i+1]-data[i])/dx;
  }
  return d;
}

////////////////////////////////////////////////////////////

v1D autocorrAbs(const v1D& t) {
  const size_t N = t.size();
  v1D ac(N, 0.); // autocorr values
  for (size_t k = 0; k<N; ++k) {
    double d = 0.;
    for (size_t i = k; i<N; ++i) {
      d += t.at(i)*t.at(i-k);
    }
    ac.at(k) = d/(double)N;
  }
  return ac;
}

////////////////////////////////////////////////////////////

v1D autocorr(const v1D& t) {
  const size_t N = t.size();
  v1D ac(N, 0.0);
  double m = average(t);
  double sum2 = 0.;
  for (size_t i = 0; i<N; ++i) {
    sum2 += (t.at(i)-m)*(t.at(i)-m);
  }
  for (size_t k = 0; k<N; ++k) {
    double sum1 = 0.;
    for (size_t i = k; i<N; ++i) {
      sum1 += (t.at(i)-m)*(t.at(i-k)-m);
    }
    ac.at(k) = sum1/sum2;
  }
  return ac;
}

////////////////////////////////////////////////////////////

CMwRandIntChooser::CMwRandIntChooser(const std::vector<double> &  weight)
: m_range( weight.size() )
{
  const size_t N = weight.size();
  if(N>(size_t)RAND_MAX)
    throw std::logic_error(std::string("CMwRandIntChooser : N > RAND_MAX"));
  // normalizing constant
  double norm = 0.0;
  for(size_t i=0; i!=N; ++i) {
    norm += weight[i];
  }
  //cerr << "[CMwRandIntChooser] norm= " << norm << endl;
  // move from double values to integers
  m_range[0]=0;
  for(size_t i=1; i!=N; ++i) {
    m_range[i] = m_range[i-1] + static_cast<int>((weight[i-1]/norm)*RAND_MAX);
    //cerr << "[CMwRandIntChooser] m_range[" << i << "]=";
    //cerr << m_range[i] << " ";
    //cerr << (double)m_range[i]/INT_MAX << endl;
  }
}


unsigned int CMwRandIntChooser::getRandInt() const
{
  const int r = rand();
  const size_t N = m_range.size();
  for(size_t i=0; i!=N-1; ++i) {
    if(r<m_range[i+1])
      return (unsigned int)i;
  }
  return (unsigned int)(N-1);
}

} // end of 'namespace mw'
