#ifndef _TRAJECTORY_GENERATOR_DDP_H_
#define _TRAJECTORY_GENERATOR_DDP_H_

#include <utils/data_type.h>
#include <utils/bezier_base.h>
// #include <global_planner/spatial_optimizer.h>
#include "mader_types.hpp"


typedef std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>> VecOfVecXd;
typedef std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd>> VecOfMatXd;

class algParam
{
    private:
        double mu;
        int maxiter;
        double tol;
        bool infeas;

        friend class ddpTrajOptimizer; // Friend Class
        friend class fwdPass; 
};

class fwdPass
{
    private:

        // high-level parameter
        int traj_order;
        int dim;
        int sys_order;
        int num_ctrlP;
        int N; // horizon

        double maxVel;
        double maxAcc;

        // defined in dynamics function
        VecOfVecXd x;  // x
        VecOfVecXd u;  // u
        VecOfVecXd c;     // c
        VecOfVecXd y; // y
        VecOfVecXd s; // s
        VecOfVecXd mu; // mu

        // defined in dynamics function
        double p;
        Eigen::VectorXd px;
        Eigen::MatrixXd pxx;
        VecOfMatXd f;     // f
        VecOfMatXd fx;
        VecOfMatXd fu;

        VecOfMatXd fxx;
        VecOfMatXd fxu;
        VecOfMatXd fuu;

        Eigen::VectorXd q;
        VecOfMatXd qx;
        VecOfMatXd qu;
        VecOfMatXd qxx;
        VecOfMatXd qxu;
        VecOfMatXd quu;

        VecOfMatXd cx;
        VecOfMatXd cu;

        // used in fx fu
        Eigen::MatrixXd barFk;
        Eigen::MatrixXd barGk;
        Eigen::MatrixXd barFkprime;
        Eigen::MatrixXd barGkprime;        
        Eigen::MatrixXd barFkpprime;
        Eigen::MatrixXd barGkpprime;   
        // used in constraints
        std::vector<double> Ek_inv;
        Eigen::VectorXd barEk_inv;

        Eigen::MatrixXd poly2bez; // berstern = poly2bez * poly;
        Eigen::MatrixXd t2tauMat; // t to tau matrix

        Eigen::MatrixXd polyt2beztau;
        Eigen::MatrixXd beztau2polyt;

        bool minvo_enabled;
        Eigen::MatrixXd polyt2minvotau;
        Eigen::MatrixXd barpoly2minvotau_v;
        Eigen::MatrixXd barpoly2minvotau_a;

        // used in cost
        Eigen::MatrixXd Pmat; 
        double w_snap;
        double Rtime;
        int time_power;
        Eigen::VectorXd x_d;  // desired terminal state
        Eigen::MatrixXd barR; // cost matrix in q 
        Eigen::MatrixXd barRprime; // cost matrix in q 
        Eigen::MatrixXd barRpprime; // cost matrix in q 


        Eigen::MatrixXd filter; //filter 
        // defined in initialroll function
        double cost;  // cost
        double costq; // without terminal
        // defined in resetfilter function
        double err;
        double logcost;
        int step;
        bool failed;
        double stepsize;
        friend class ddpTrajOptimizer; // Friend Class
        Eigen::Matrix<double, 6, 6> polyt2minvotau6_;
        Eigen::Matrix<double, 5, 6> polyt2minvotau_v6_;
        Eigen::Matrix<double, 4, 6> polyt2minvotau_a6_;

        Eigen::VectorXd jerkCost;
        Eigen::MatrixXd comMat;

        double reg_exp_base;
    public:

        // for call in forward pass or in the function of class fwdPass
        fwdPass(){            
            comMat = Eigen::MatrixXd::Zero(18,18);
            comMat << 1.0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                       0,   0,   0,   0,   0,   0, 1.0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                       0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 1.0,   0,   0,   0,   0,   0,
                       0, 1.0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                       0,   0,   0,   0,   0,   0,   0, 1.0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                       0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 1.0,   0,   0,   0,   0,
                       0,   0, 1.0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                       0,   0,   0,   0,   0,   0,   0,   0, 1.0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                       0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 1.0,   0,   0,   0,
                       0,   0,   0, 1.0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                       0,   0,   0,   0,   0,   0,   0,   0,   0, 1.0,   0,   0,   0,   0,   0,   0,   0,   0,
                       0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 1.0,   0,   0,
                       0,   0,   0,   0, 1.0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                       0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 1.0,   0,   0,   0,   0,   0,   0,   0,
                       0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 1.0,   0,
                       0,   0,   0,   0,   0, 1.0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                       0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 1.0,   0,   0,   0,   0,   0,   0,
                       0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 1.0;
}
        ~fwdPass(){}
        void setPmat();
        void setPmat(Eigen::MatrixXd Pmat_input);
        void time2barFkbarGk(double time);
        void time2barFkprimebarGkprime(double time);
        void time2barR(double time);
        void t2tau(double time);
        Eigen::VectorXd  computenextx(Eigen::VectorXd xnew, Eigen::VectorXd unew);
        Eigen::VectorXd  computec(Eigen::VectorXd xnew, Eigen::VectorXd unew, int time_id, const decomp_cvx_space::FlightCorridor &corridor);
        Eigen::VectorXd  computecminvo(const Eigen::VectorXd& xnew, const Eigen::VectorXd& unew, 
                                        int time_id, const decomp_cvx_space::FlightCorridor &corridor);
        // Eigen::VectorXd  computecminvomixed(const Eigen::VectorXd& xnew, const Eigen::VectorXd& unew, 
                                        // int time_id, const decomp_cvx_space::FlightCorridor &corridor);

        double  computep(Eigen::VectorXd& xnew, Eigen::VectorXd& xnew_d);
        double  computeq(Eigen::VectorXd& xnew, Eigen::VectorXd& unew);


        // values and derivatives along entire traj
        void computeall(const decomp_cvx_space::FlightCorridor &corridor);
        void computeprelated();
        void computefrelated();
        void computeqrelated();
        void computecrelated(const decomp_cvx_space::FlightCorridor &corridor);
        void computecrelatedminvo( const decomp_cvx_space::FlightCorridor &corridor);
        // void computecrelatedminvomixed( const decomp_cvx_space::FlightCorridor &corridor);

        // for call in optimizer
        void initialroll(const decomp_cvx_space::FlightCorridor  &corridor);
        void resetfilter( algParam alg );
        void finalroll();

        // for call in backward pass;
        VecOfVecXd getx() {return x;}; 
        VecOfVecXd getu() {return u;};
        VecOfVecXd getc() {return c;};

        VecOfVecXd gety() {return y;};
        VecOfVecXd gets() {return s;};
        VecOfVecXd getmu() {return mu;}; 

        double getp(){return p;};
        Eigen::VectorXd getpx() {return px;};
        Eigen::MatrixXd getpxx() {return pxx;};

        VecOfMatXd getfx() {return fx;};
        VecOfMatXd getfu() {return fu;};
        VecOfMatXd getfxx() {return fxx;};
        VecOfMatXd getfxu() {return fxu;};
        VecOfMatXd getfuu() {return fuu;};

        VecOfMatXd getqx() {return qx;};
        VecOfMatXd getqu() {return qu;};
        VecOfMatXd getqxx() {return qxx;};
        VecOfMatXd getqxu() {return qxu;};
        VecOfMatXd getquu() {return quu;};

        VecOfMatXd getcx() {return cx;};
        VecOfMatXd getcu() {return cu;};

        void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);

};

class bwdPass
{
    private:
        // defined in dynamic function
        VecOfVecXd ky; // ky
        VecOfMatXd Ky; // Ky
        VecOfVecXd ks; // ks
        VecOfMatXd Ks; // Ks
        VecOfVecXd ku; // ku
        VecOfMatXd Ku; // Ku

        // defined in resetreg function
        double reg;
        bool failed;
        double recovery;

        double opterr;
        Eigen::Vector2d dV;
        friend class ddpTrajOptimizer; // Friend Class

    public:
        void resetreg();
        void initreg(double regvalue);

};




class ddpTrajOptimizer 
{
    private:
        double ddpobj;
        // double snapobj;
        Eigen::MatrixXd PolyCoeff;
        Eigen::MatrixXd BezCoeff;

        Eigen::VectorXd PolyTime;
        double compTime;
        VecOfVecXd PolyTimeTraj;
        VecOfMatXd PolyCoeffTraj;
        std::vector<double> costTraj;
        std::vector<double> costqTraj;
        // include this in class, so that no need to pass
        algParam alg;
        fwdPass fp;
        bwdPass bp;

        int iter_used;
        Eigen::VectorXd jerkCost;
        mt::PieceWisePol pwp_out_;

    public:
        ddpTrajOptimizer(){}
        ~ddpTrajOptimizer(){}

        int polyCurveGeneration( 
        const decomp_cvx_space::FlightCorridor &corridor,
        // const Eigen::MatrixXd &MQM_u,
        // const Eigen::MatrixXd &MQM_l,
        const Eigen::MatrixXd &pos,
        const Eigen::MatrixXd &vel,
        const Eigen::MatrixXd &acc,
        const Eigen::MatrixXd &jer,
        const double max_vel, 
        const double max_acc,
        const double max_jer,
        const mt::PieceWisePol &pwp_init,
        const double w_snap,
        const double w_terminal,
        const double w_time,
        const int iter_max,
        bool &infeas,
        bool zero_init_flag,
        bool line_init_flag,
        bool &line_failed,
        int time_power,
        bool minvo_flag);

        void backwardpass(const decomp_cvx_space::FlightCorridor &corridor);
        void forwardpass( const decomp_cvx_space::FlightCorridor &corridor);        

        Eigen::MatrixXd poly2bezFunc();        
        Eigen::MatrixXd bez2polyFunc();        

        void sysparam2polyFunc();        
        
        Eigen::MatrixXd getPolyCoeff()
        {
            return PolyCoeff;
        };

        Eigen::MatrixXd getBezCoeff()
        {
            return BezCoeff;
        };        

        Eigen::VectorXd getPolyTime()
        {
            return PolyTime;
        };

        double getDDPObjective()
        {
            return ddpobj;
        };

        
        // double getSnapObjective()
        // {
        //     return snapobj;
        // };


        double getCompTime()
        {
            return compTime;
        };   

        double getTerminalNorm()
        {
            return (fp.x.back() - fp.x_d).transpose() * (fp.x.back() - fp.x_d);
        }     

        int getIterUsed()
        {
            return iter_used;
        }

        double getJerkCost()
        {
            return jerkCost.sum();
        }
        void generatePwpOut(mt::PieceWisePol& pwp_out, 
                            std::vector<mt::state>& traj_out, double t_start, double dc);
};

struct resultsStruct
{
    int num_of_poly;
    double compTime; // from optimizer
    double allocTime;
    double initallocTime;
};

#endif