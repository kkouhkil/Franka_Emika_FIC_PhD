#include "../include/CartesianControl.hpp"
#include "../include/DynamicsControl.hpp"
#include "../include/JointData.hpp"
#include <CosimaUtilities/SkewMatrix.hpp>
#include <CosimaUtilities/TorqueSaturation.hpp>

class CartImpedanceControl : public RTT::TaskContext {
public:
    CartImpedanceControl(std::string const &name);

    bool configureHook();
    bool startHook();
    void updateHook();
    void stopHook();
    void cleanupHook();

    void compute();
    void setGains(float gainP, float gainD, float orientation_gainP, float orientation_gainD);
	void setGains_anisotropic(const Eigen::VectorXf& Stiffness_Des, const Eigen::VectorXf& Damping_Des);

    void addChain(int dof);
    void preparePorts();
    void displayCurrentState();
    double getSimulationTime();

// // // // // // // // // // // // // // // //// // // // // NEW PART OF THE CODE // // // // // // // // // // // // // // // //// // // // // - 11.02.2019
    void quatToEuler(const double qW, const double qX, const double qY, const double qZ, Eigen::VectorXf &EulerAngle);
    void calculate_K_Max_Beta_Pos_Ori(Eigen::VectorXf &F_MAX, Eigen::VectorXf &Max_Disp, Eigen::VectorXf &Tau_MAX, Eigen::VectorXf &Max_Disp_Ori, Eigen::MatrixXf &K_Max_System, Eigen::VectorXf &Beta ,Eigen::MatrixXf &K_Max_System_Ori, Eigen::VectorXf &Beta_Ori, int i);
    void get_K_Total_Pos_Ori(Eigen::VectorXf &Beta, Eigen::MatrixXf &K_Var_Pos, Eigen::MatrixXf &K_Total_Des_Pos, Eigen::VectorXf &Beta_Ori, Eigen::MatrixXf &K_Var_Ori, Eigen::MatrixXf &K_Total_Des_Ori, int i);
    void calculate_F_Max (Eigen::VectorXf &Max_Disp, Eigen::VectorXf &F_Max, int i);
    void calculate_Tau_Max(Eigen::VectorXf &Max_Disp_Ori, Eigen::VectorXf &Tau_Max, int i);
    void boxShaping_Online(Eigen::VectorXf &Max_Disp, int i);
    void calculation_K1_K2_Kd_Max_PMax_PMid_Max_Dist_Vec(Eigen::MatrixXf &K1, Eigen::MatrixXf &K2, Eigen::MatrixXf &Kd_Max, Eigen::VectorXf &Pmax, Eigen::VectorXf &Pmid, Eigen::VectorXf &Max_Dist_Vec, Eigen::VectorXf &Beta, Eigen::MatrixXf &K_Max_System, Eigen::VectorXf &Max_Disp, int i);
    void calculation_K1_K2_Kd_Max_PMax_PMid_Max_Dist_Vec_Ori(Eigen::MatrixXf &K1_Ori, Eigen::MatrixXf &K2_Ori, Eigen::MatrixXf &Kd_Max_Ori, Eigen::VectorXf &Pmax_Ori, Eigen::VectorXf &Pmid_Ori, Eigen::VectorXf &Max_Dist_Vec_Ori, Eigen::VectorXf &Beta_Ori, Eigen::MatrixXf &K_Max_System_Ori, int i);
    void calculation_Kc_DampConstLiner_FMidCtrl_FkCtrl(Eigen::MatrixXf &Kc, Eigen::MatrixXf &DampConstLinear, Eigen::VectorXf &FMidCtrl, Eigen::VectorXf &FkCtrl, int i);
    void calculation_Kc_DampConstAngular_FMidCtrl_FkCtrl_Ori(Eigen::MatrixXf &Kc_Ori, Eigen::MatrixXf &DampConstAngular, Eigen::VectorXf &FMidCtrl_Ori, Eigen::VectorXf &FkCtrl_Ori, int i);
    void trajectory_execution(Eigen::VectorXf &DesPos, Eigen::VectorXf &EndEff_VelLinear, int i);

    template <typename T> static int sgn(const T &val, const T &eps = T(0)) {
            return (T(eps) < val) - (val < T(-eps));
        }

    // // // // // // // // // // // // // // // //// // // // // - 11.02.2019 // // // // // // // // // // // // // // // //// // // // //

private:	
	
    std::unique_ptr<CartesianControl> cc;
    std::unique_ptr<DynamicsControl> dc;
    std::unique_ptr<JointData> jd;

    std::unique_ptr<CosimaUtilities::Pseudoinverse<Eigen::MatrixXf> > pinv;
    std::unique_ptr<CosimaUtilities::Pseudoinverse<Eigen::MatrixXf> > nullSpace;
    std::unique_ptr<CosimaUtilities::TorqueSaturation<Eigen::VectorXf> > trqSat;

    float thresh = 1e-4;

    int numOfEndEffectors;

	// TEXT FILE REQUIRED INPUTS
    int numTrajAct, numTraj_CircTraj_Act , numBoxShapingAct;
    double numTrajSize, numTrajOmega, numCircTrajSize;
    double numVirtPosBound_x, numVirtPosBound_y, numVirtPosBound_z, numVirtOriBound_x, numVirtOriBound_y, numVirtOriBound_z, numDesPos_x, numDesPos_y, numDesPos_z;
    //
	
	double gainP, gainD, orientation_gainP, orientation_gainD;
	double currentTime, timeInitial,time0, time1, time2, deltaTime;
	double omegaTraj, lineDispTrajCoef, circRadiusTraj;
    double slopeFmax, slopeTauMax_Ori, onlineBoxShapingSwitch;
	double trajStartTime;
	double trajectory_activation;

	double potentialEnergy, kineticEnergy;

	Eigen::VectorXf total_system_energy, energy_vector;
	Eigen::MatrixXf lambda_c;
	Eigen::VectorXf Stiffness_Des, Damping_Des;
	Eigen::VectorXf Fd_Pos, Fd_Ori;


	Eigen::VectorXf curPos3x1, desPos3x1, curOri3x1, desOri3x1, curVel3x1, desVel3x1;
 	Eigen::VectorXf position_error, angular_error;

	Eigen::VectorXf posError_oriError_6x1;
	Eigen::VectorXf velLinearError, velAngularError, endEffectorVelError;

	Eigen::MatrixXf kConstPos, kVarPos, kTotalDesPos;
    Eigen::MatrixXf kConstOri, kVarOri, kTotalDesOri;
    
	Eigen::MatrixXf dConstLinear, dVarPos;
    Eigen::MatrixXf dConstAngular, dVarOri;

	Eigen::MatrixXf kTotal6x6;
  	Eigen::MatrixXf dTotal6x6, dVar6x6, dConst6x6, dCritical6x6;  
    
	Eigen::VectorXf curJntPos7x1, desJntPos7x1;
    Eigen::VectorXf curJntVel7x1, desJntVel7x1;

   	Eigen::MatrixXf kConstNull, kVarNull;
    Eigen::MatrixXf dConstNull, dVarNull;

    Eigen::VectorXf preEndEff_PosError3x1;
    Eigen::VectorXf preEndEff_VelError3x1;

    Eigen::VectorXf preEndEff_OriError3x1;
    Eigen::VectorXf preEndEff_VelOriError3x1;

    Eigen::VectorXf nullSpaceCtrl;
    Eigen::VectorXf torqueSaturation;

    Eigen::MatrixXf KdMax;
    Eigen::MatrixXf KdMax_Ori;

	Eigen::VectorXf maxVel3x1;
    Eigen::VectorXf Pmax3x1;
    Eigen::VectorXf Pmid3x1;

    Eigen::MatrixXf K1_3x3;
    Eigen::MatrixXf K2_3x3;
    Eigen::MatrixXf Kc_3x3;
    Eigen::VectorXf Fk;
    Eigen::VectorXf F_midPointErrorCtrl;
    Eigen::VectorXf maxDistVec_3x1;

    Eigen::VectorXf maxVel3x1_Ori;
    Eigen::VectorXf Pmax3x1_Ori;
    Eigen::VectorXf Pmid3x1_Ori;
    Eigen::MatrixXf K1_Ori_3x3;
    Eigen::MatrixXf K2_Ori_3x3;
    Eigen::MatrixXf Kc_Ori_3x3;
    Eigen::VectorXf Fk_Ori;
    Eigen::VectorXf F_midPointErrorCtrl_Ori;
    Eigen::VectorXf maxDistVec_Ori_3x1;

    // PROPOSED METHOD 6X1 MAXIMUM DISPLACEMENT VECTOR FOR ENERGY CALCULATION
    Eigen::VectorXf maxDistVec_6x1;

    Eigen::VectorXf Fmax;
    Eigen::VectorXf maxDisp;
    Eigen::VectorXf beta;
    Eigen::MatrixXf KmaxSystem;

    Eigen::VectorXf TauMax;
    Eigen::VectorXf maxDisp_Ori;
    Eigen::VectorXf beta_Ori;
    Eigen::MatrixXf KmaxSystem_Ori;

	Eigen::VectorXf desPos_curPos;
	Eigen::VectorXf Fctrl_XYZ;
	//Eigen::VectorXf endEff_force_vec;


 	double sinr_cosp, cosr_cosp, roll, sinp, pitch, siny_cosp, cosy_cosp, yaw;

	// NEW PART OF THE CODE // 30.01.2019
	Eigen::VectorXf Fdamping;
	/////////////////////////////////////

	// Link 3 end effector code
    //RTT::InputPort<Eigen::MatrixXf> in_jacobianLink3_port;
    //RTT::FlowStatus in_jacobianLink3_Flow;
    //Eigen::MatrixXf JacLink3;
    //RTT::InputPort<Eigen::VectorXf> in_DirectKLink3_port;
    //RTT::FlowStatus in_DirectKLink3_Flow;
    //Eigen::VectorXf DirectKLink3;

	Eigen::VectorXf estimated_torque_var;

    RTT::OutputPort<Eigen::VectorXf> out_estimated_torque_port;
    RTT::OutputPort<Eigen::VectorXf> out_variable_stiffness_pos_ori_port;
    RTT::OutputPort<Eigen::VectorXf> out_posError_oriError_port;
    RTT::OutputPort<Eigen::VectorXf> out_energyTank_port;
    RTT::OutputPort<Eigen::VectorXf> out_critical_damping_pos_ori_port;
    RTT::OutputPort<Eigen::VectorXf> out_posError_desTraj_curPos_port;
    RTT::OutputPort<Eigen::VectorXf> out_varDamping_port;
    RTT::OutputPort<Eigen::VectorXf> out_inOutPower_dampPowerCoef_port;
    RTT::OutputPort<Eigen::VectorXf> out_endEffectorVelError_port;
	RTT::OutputPort<Eigen::VectorXf> out_systemEnergy_port;

	RTT::InputPort<Eigen::MatrixXf> in_jacobian_port;

	
	// NEW PART OF THE CODE // 28.02.2019
	RTT::OutputPort<Eigen::VectorXf> out_desPos_curPos_port;
	RTT::OutputPort<Eigen::VectorXf> out_Fctrl_port;
	//RTT::OutputPort<Eigen::VectorXf> out_endEff_force_port;
    
	//// //// //// //// //// //// BUNCH OF JUNKS! //// //// //// //// //// ////
	double x_tDot, x_t, energyUpperBound, energyLowerBound, T, sigma, wDot, x_t_New, x_tDot_New, T_New, inOutPower, dampPowerCoef, PcMax, Pc, D0, cSlopeFmax;
	Eigen::MatrixXf dCriticalPos;
    Eigen::MatrixXf dCriticalOri;
	Eigen::VectorXf inOutPowerVec;
	Eigen::VectorXf desTraj3x1;
    Eigen::VectorXf desTrajVel3x1;
    Eigen::VectorXf posError_desTraj_curPos_6x1;
    Eigen::VectorXf inOutPower_dampPowerCoef3x1;
	Eigen::VectorXf omega, omega_New;
    Eigen::VectorXf phi, phi_New;
    Eigen::MatrixXf Ktraj_3x3;
	Eigen::MatrixXf MK6x6;
    Eigen::MatrixXf eigenVector6x6;
    Eigen::MatrixXf eigenValue6x6;
    Eigen::MatrixXf eigenVaule6x6Sqrt;
	Eigen::VectorXf psi, psi_New;
	Eigen::VectorXf varStiffVec;
    Eigen::VectorXf energyTankVec;
    Eigen::VectorXf dCriticalPosOri;
	Eigen::VectorXf estimated_Fext;
    Eigen::VectorXf momentumVec;
};
