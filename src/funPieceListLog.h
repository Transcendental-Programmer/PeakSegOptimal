/* -*- compile-command: "R CMD INSTALL .." -*- */

#define ERROR_MIN_MAX_SAME 1
#include <list>

// NOTE: please only define prototypes in this file (do not define
// methods directly -- instead define them in funPieceList.cpp). This
// is more compatible with R's Makefile, which automatically
// recompiles object files when there are changes to *.cpp but not *.h
// files.

// Add new base class for loss pieces.
class LossPieceBase {
 public:
    virtual double getCost(double x) = 0;
    virtual double getDeriv(double x) = 0;
    virtual double argmin() = 0;
    virtual void print() = 0;
    virtual ~LossPieceBase() {}
};

// Modify PoissonLossPieceLog to inherit from LossPieceBase.
class PoissonLossPieceLog : public LossPieceBase {
 public:
  double Linear;
  double Log;
  double Constant;
  double min_log_mean;
  double max_log_mean;
  int data_i;
  double prev_log_mean;
  PoissonLossPieceLog();
  PoissonLossPieceLog(double li, double lo, double co, double m, double M, int i, double);
  double argmin();
  double argmin_mean();
  void print();
  double get_smaller_root(double);
  double get_larger_root(double);
  bool has_two_roots(double);
  double getCost(double mean);
  double getDeriv(double);
  double PoissonLoss(double);
  double PoissonDeriv(double);
};

// Add new NormalLossPiece class for normal loss
class NormalLossPiece : public LossPieceBase {
 public:
  double Quadratic;
  double Linear;
  double Constant;
  double min_mean;
  double max_mean;
  int data_i;
  double prev_mean;
  NormalLossPiece();
  NormalLossPiece(double q, double l, double c, double min_m, double max_m, int i, double prev);
  double argmin();
  void print();
  double get_smaller_root(double);
  double get_larger_root(double);
  bool has_two_roots(double);
  double getCost(double mean);
  double getDeriv(double);
};

typedef std::list<PoissonLossPieceLog> PoissonLossPieceListLog;
typedef std::list<NormalLossPiece> NormalLossPieceList;

class PiecewisePoissonLossLog {
 public:
  PoissonLossPieceListLog piece_list;
  void set_to_min_less_of(PiecewisePoissonLossLog *, int);
  void set_to_min_more_of(PiecewisePoissonLossLog *, int);
  void set_to_min_env_of
    (PiecewisePoissonLossLog *, PiecewisePoissonLossLog *, int);
  int check_min_of(PiecewisePoissonLossLog *, PiecewisePoissonLossLog *);
  void push_min_pieces
    (PiecewisePoissonLossLog *, PiecewisePoissonLossLog *,
     PoissonLossPieceListLog::iterator, PoissonLossPieceListLog::iterator, int);
  void push_piece(PoissonLossPieceListLog::iterator, double, double);
  void add(double Linear, double Log, double Constant);
  void multiply(double);
  void print();
  void set_prev_seg_end(int prev_seg_end);
  void findMean(double mean, int *seg_end, double *prev_log_mean);
  double findCost(double mean);
  void Minimize
    (double *best_cost,
     double *best_mean,
     int *data_i,
     double *prev_log_mean);
  
  // New method for unconstrained minimization:
  void set_to_unconstrained_min_of(PiecewisePoissonLossLog *input, int verbose);
};

// Add new class for handling normal loss functions
class PiecewiseNormalLoss {
 public:
  NormalLossPieceList piece_list;
  void set_to_min_less_of(PiecewiseNormalLoss *, int);
  int check_min_of(PiecewiseNormalLoss *, PiecewiseNormalLoss *);
  void push_min_pieces
    (PiecewiseNormalLoss *, PiecewiseNormalLoss *,
     NormalLossPieceList::iterator, NormalLossPieceList::iterator, int);
  void push_piece(NormalLossPieceList::iterator, double, double);
  void print();
  void add(double Quadratic, double Linear, double Constant);
  void multiply(double);
  void set_prev_seg_end(int prev_seg_end);
  void findMean(double mean, int *seg_end, double *prev_mean);
  double findCost(double mean);
  void Minimize
    (double *best_cost,
     double *best_mean,
     int *data_i,
     double *prev_mean);
};

bool sameFuns(PoissonLossPieceListLog::iterator, PoissonLossPieceListLog::iterator);
bool sameFunsNormal(NormalLossPieceList::iterator, NormalLossPieceList::iterator);
