#include "../include/CartImpedanceControl.hpp"
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>

#include <iostream>
#include <fstream>
#include <string>

CartImpedanceControl::CartImpedanceControl(std::string const &name) : RTT::TaskContext(name) {
    
	//OPERATIONS PREPARATION
    addOperation("setGains", &CartImpedanceControl::setGains, this).doc("set Gains");
	addOperation("setGains_anisotropic", &CartImpedanceControl::setGains_anisotropic, this).doc("set Gains anisotropic");
    addOperation("addChain", &CartImpedanceControl::addChain, this).doc("add Chain");
    addOperation("preparePorts", &CartImpedanceControl::preparePorts, this).doc("preparePorts");
    addOperation("displayCurrentState", &CartImpedanceControl::displayCurrentState, this).doc("display current state");
    addProperty("thresh", thresh);
    

	// DEFINING THE USE OF CARTESIAN CTRL, DYNAMICS CTRL AND JOINT DATA
	cc = std::make_unique<CartesianControl>(this);
    dc = std::make_unique<DynamicsControl>(this);
    jd = std::make_unique<JointData>(this);

	// NUMBER OF END-EFFECTOR(S) INITIALIZATION
	numOfEndEffectors = 0;

    //INITIALIZING AND READING ELEMENTS FROM TXT FILE
    numTrajAct = 0; numTraj_CircTraj_Act = 0; numBoxShapingAct = 0; numTrajSize = 0; numTrajOmega = 0; numCircTrajSize = 0;
    numVirtPosBound_x = 0; numVirtPosBound_y = 0; numVirtPosBound_z = 0;
    numVirtOriBound_x = 0; numVirtOriBound_y = 0; numVirtOriBound_z = 0;
    numDesPos_x = 0; numDesPos_y = 0; numDesPos_z = 0;

    std::fstream inFile;
    inFile.open("/home/joshua/workspaces/keyhan_ws/src/cartesian_var_imp_ctrl/numbers.txt");
    inFile >> numTrajAct >> numTraj_CircTraj_Act >> numBoxShapingAct >>  numTrajSize >> numTrajOmega >> numCircTrajSize >>
                            numVirtPosBound_x >> numVirtPosBound_y >> numVirtPosBound_z >> numVirtOriBound_x >> numVirtOriBound_y >> numVirtOriBound_z >>
                            numDesPos_x >> numDesPos_y >> numDesPos_z;

	// SETTING CONSTANT GAINS - ELEMENTS
    gainP = 0; gainD = 0; D0 = 0;
    orientation_gainP = 0; orientation_gainD = 0;
    //wDot = 0;

	// CALCULATION OF DELTA-TIME - ELEMENTS
    timeInitial = 0;time1 = 0; time2 = 0;
    deltaTime = 0;

	// TASK-SPACE - ELEMENTS
    curPos3x1 = Eigen::VectorXf (3); curPos3x1.setZero();
    desPos3x1 = Eigen::VectorXf (3); desPos3x1.setZero() ; //desPos3x1 = cc->des_poses_var.head<3>();
    curOri3x1 = Eigen::VectorXf (3); curOri3x1.setZero();
    desOri3x1 = Eigen::VectorXf (3); desOri3x1.setZero();
    curVel3x1 = Eigen::VectorXf (3); curVel3x1.setZero();
    desVel3x1 = Eigen::VectorXf (3); desVel3x1.setZero();

    velLinearError = Eigen::VectorXf (3); velLinearError.setZero();
	endEffectorVelError = Eigen::VectorXf (6); endEffectorVelError.setZero();
    velAngularError = Eigen::VectorXf (3); velAngularError.setZero();

    position_error = Eigen::VectorXf (3); position_error.setZero();
    angular_error = Eigen::VectorXf (3); angular_error.setZero();

	posError_oriError_6x1 = Eigen::VectorXf (6); posError_oriError_6x1.setZero();

    kConstPos = Eigen::MatrixXf(3, 3); kConstPos.setIdentity();
    addProperty("kConstPos_yy", kConstPos(1,1));
    kVarPos = Eigen::MatrixXf(3, 3); kVarPos.setIdentity();
	kTotalDesPos = Eigen::MatrixXf(3, 3); kTotalDesPos.setIdentity();
    kConstOri = Eigen::MatrixXf(3, 3); kConstOri.setIdentity();
    kVarOri = Eigen::MatrixXf(3, 3); kVarOri.setIdentity();
	kTotalDesOri = Eigen::MatrixXf(3, 3); kTotalDesOri.setIdentity();
	varStiffVec = Eigen::VectorXf (6); varStiffVec.setZero();

    dConstLinear = Eigen::MatrixXf (3,3); dConstLinear.setIdentity();
    dVarPos = Eigen::MatrixXf(3, 3); dVarPos.setIdentity();
    dConstAngular = Eigen::MatrixXf (3, 3); dConstAngular.setIdentity();
    dVarOri = Eigen::MatrixXf(3, 3); dVarOri.setIdentity();

    kTotal6x6 = Eigen::MatrixXf (6,6); kTotal6x6.setIdentity();
    dTotal6x6 = Eigen::MatrixXf (6,6); dTotal6x6.setIdentity();
	dConst6x6 = Eigen::MatrixXf (6,6); dConst6x6.setIdentity();
	dCritical6x6 = Eigen::MatrixXf (6,6); dCritical6x6.setIdentity();

    // DAMPING FORCES FOR LINEAR AND ANGULAR VELOCITIES - REQUIRED ELEMENTS
    Fd_Pos = Eigen::VectorXf(3); Fd_Pos.setZero();
    Fd_Ori = Eigen::VectorXf(3); Fd_Ori.setZero();

	// NULL-SPACE - ELEMENTS
    curJntPos7x1 = Eigen::VectorXf (7); curJntPos7x1.setZero();
    desJntPos7x1 = Eigen::VectorXf (7); desJntPos7x1.setZero();
    curJntVel7x1 = Eigen::VectorXf (7); curJntVel7x1.setZero();
    desJntVel7x1 = Eigen::VectorXf (7); desJntVel7x1.setZero();

    // NULL-SPACE - STIFFNESS - ELEMENTS
    kConstNull = Eigen::MatrixXf (7,7); kConstNull.setIdentity();
    kVarNull = Eigen::MatrixXf (7,7); kVarNull.setIdentity();

    // NULL-SPACE - DAMPING - ELEMENTS
    dConstNull = Eigen::MatrixXf (7,7); dConstNull.setIdentity();
    dVarNull = Eigen::MatrixXf (7,7); dVarNull.setIdentity();

    // CALCULATION OF PREVIOUS VALUES - ELEMENTS
    preEndEff_PosError3x1 = Eigen::VectorXf(3); preEndEff_PosError3x1.setZero();
    preEndEff_OriError3x1 = Eigen::VectorXf(3); preEndEff_OriError3x1.setZero();
    preEndEff_VelError3x1 = Eigen::VectorXf(3); preEndEff_VelError3x1.setZero();
    preEndEff_VelOriError3x1 = Eigen::VectorXf(3); preEndEff_VelOriError3x1.setZero();

    // NULL-SPACE CTRL PLUS TORQUE SATURATION - ELEMENTS
    nullSpaceCtrl = Eigen::VectorXf (7); nullSpaceCtrl.setZero();
    torqueSaturation = Eigen::VectorXf (7); torqueSaturation.setZero();

    // KdMAX - ELEMENTS
    KdMax = Eigen::MatrixXf(3, 3); KdMax.setIdentity();
    KdMax_Ori = Eigen::MatrixXf(3, 3); KdMax_Ori.setIdentity();

    // SYSTEM ENERGY - POTENTIAL AND KINETIC INITIALIZATION
    potentialEnergy = 0; kineticEnergy = 0;
	total_system_energy = Eigen::VectorXf(3); total_system_energy.setZero();
	energy_vector = Eigen::VectorXf(3); energy_vector.setZero();

    // PROPOSED METHOD - POSITION
    Pmax3x1 = Eigen::VectorXf(3); Pmax3x1.setZero();
    Pmid3x1 = Eigen::VectorXf(3); Pmid3x1.setZero();
    maxVel3x1 = Eigen::VectorXf(3); maxVel3x1.setZero();
    maxDistVec_3x1 = Eigen::VectorXf(3); maxDistVec_3x1.setZero();
    K1_3x3 = Eigen::MatrixXf(3,3); K1_3x3.setIdentity();
    K2_3x3 = Eigen::MatrixXf(3,3); K2_3x3.setIdentity();
    Kc_3x3 = Eigen::MatrixXf(3,3); Kc_3x3.setZero();
    F_midPointErrorCtrl = Eigen::VectorXf(3); F_midPointErrorCtrl.setZero();
    Fk = Eigen::VectorXf(3); Fk.setZero();

    // PROPOSED METHOD - ORIENTATION - ELEMENTS
    maxVel3x1_Ori = Eigen::VectorXf(3); maxVel3x1_Ori.setZero();
    Pmax3x1_Ori = Eigen::VectorXf(3);Pmax3x1_Ori.setZero();
    Pmid3x1_Ori = Eigen::VectorXf(3);Pmid3x1_Ori.setZero();
    maxDistVec_Ori_3x1 = Eigen::VectorXf(3); maxDistVec_Ori_3x1.setZero();
    K1_Ori_3x3 = Eigen::MatrixXf(3,3); K1_Ori_3x3.setIdentity();
    K2_Ori_3x3 = Eigen::MatrixXf(3,3); K2_Ori_3x3.setIdentity();
    Kc_Ori_3x3 = Eigen::MatrixXf(3,3); Kc_Ori_3x3.setIdentity();
    F_midPointErrorCtrl_Ori = Eigen::VectorXf(3); F_midPointErrorCtrl_Ori.setZero();
    Fk_Ori = Eigen::VectorXf(3); Fk_Ori.setZero();

	// NEW PART OF THE CODE
    Fmax = Eigen::VectorXf(3); Fmax.setZero();

    maxDisp = Eigen::VectorXf(3); //maxDisp.setZero();
    beta = Eigen::VectorXf(3); beta.setZero();
    addProperty("beta_y", beta[1]);
    KmaxSystem = Eigen::MatrixXf(3,3); KmaxSystem.setIdentity();
     addProperty("KMax_yy", KmaxSystem(1,1));
	//maxDisp.setConstant(0.2);
    maxDisp[0] = numVirtPosBound_x;
    maxDisp[1] = numVirtPosBound_y;
    maxDisp[2] = numVirtPosBound_z;


 // NEW PART OF THE CODE - Fmax_Ori
    TauMax = Eigen::VectorXf(3); TauMax.setZero();
    maxDisp_Ori = Eigen::VectorXf(3); //maxDisp_Ori.setZero();
    beta_Ori = Eigen::VectorXf(3); beta_Ori.setZero();
    KmaxSystem_Ori = Eigen::MatrixXf(3,3); KmaxSystem_Ori.setIdentity();
    //maxDisp_Ori.setConstant(0.2);
    maxDisp_Ori[0] = numVirtOriBound_x;
    maxDisp_Ori[1] = numVirtOriBound_y;
    maxDisp_Ori[2] = numVirtOriBound_z;

    // Fmax SLOPE INITIALIZATION
    slopeFmax = 0; slopeTauMax_Ori = 0;

    // ONLINE BOX SHAPING SWITCH
    onlineBoxShapingSwitch = numBoxShapingAct;
    time0 = 0;

    // STIFNESS AND DAMPING VECTORS FOR ANISOTROPIC VALUE SETS
    Stiffness_Des = Eigen::VectorXf(6); Stiffness_Des.setZero();
    Damping_Des = Eigen::VectorXf(6); Damping_Des.setZero();

	// TRAJECTORY EXECUTION - ELEMENTS
    trajectory_activation = numTrajAct;
	lineDispTrajCoef = numTrajSize; //2.50; //trajectory motion in cm from desired point // 25 * ((0.1/2)/10)
	circRadiusTraj = numCircTrajSize;
	omegaTraj = numTrajOmega;
	trajStartTime = 5;

	desPos_curPos = Eigen::VectorXf(6); desPos_curPos.setZero();
	Fctrl_XYZ = Eigen::VectorXf(6); Fctrl_XYZ.setZero();

	sinr_cosp = 0; cosr_cosp = 0; roll = 0; sinp = 0; pitch = 0; siny_cosp = 0; cosy_cosp = 0; yaw = 0;

	Fdamping = Eigen::VectorXf(3); Fdamping.setZero();

////////////////////////////////////////////////////////////////////////

	// DO NOT NEED TEHM, BUT CANNOT BE DELETED!!!
    inOutPowerVec = Eigen::VectorXf (6); inOutPowerVec.setZero();
    momentumVec = Eigen::VectorXf (6); momentumVec.setZero();
	Ktraj_3x3 = Eigen::MatrixXf(3,3); Ktraj_3x3.setIdentity(); 

}

void CartImpedanceControl::addChain(int dof) {
    numOfEndEffectors++;

    cc->addChain(dof);
    dc->addChain(dof);
    jd->addChain(dof);
    estimated_torque_var.resize(dof + estimated_torque_var.size());
    estimated_torque_var.setZero();
}
void CartImpedanceControl::preparePorts() {
    cc->preparePorts();
    dc->preparePorts();
    jd->preparePorts();

    Eigen::MatrixXf tmp2(6 * numOfEndEffectors, 6 * numOfEndEffectors);
    pinv = std::make_unique<CosimaUtilities::Pseudoinverse<Eigen::MatrixXf> >(tmp2);

    Eigen::MatrixXf tmp3(7 * numOfEndEffectors, 7 * numOfEndEffectors);
    nullSpace = std::make_unique<CosimaUtilities::Pseudoinverse<Eigen::MatrixXf> >(tmp3);

    Eigen::VectorXf tmp4(7 * numOfEndEffectors);
    tmp4.setZero();
    trqSat = std::make_unique<CosimaUtilities::TorqueSaturation<Eigen::VectorXf> > (0.4f, tmp4);
 
    out_estimated_torque_port.setName("out_estimated_torque_port");
    out_estimated_torque_port.doc("Output port for sending estimated force values");
    out_estimated_torque_port.setDataSample(estimated_torque_var);
    ports()->addPort(out_estimated_torque_port);

    out_variable_stiffness_pos_ori_port.setName("out_variable_stiffness_pos_ori_port");
    out_variable_stiffness_pos_ori_port.doc("Output port for printing variable stiffness Pos and Ori");
    out_variable_stiffness_pos_ori_port.setDataSample(varStiffVec);
    ports()->addPort(out_variable_stiffness_pos_ori_port);

    out_posError_oriError_port.setName("out_posError_oriError_port");
    out_posError_oriError_port.doc("Output port for printing Pos Error and Ori Error");
    out_posError_oriError_port.setDataSample(posError_oriError_6x1);
    ports()->addPort(out_posError_oriError_port);

    out_critical_damping_pos_ori_port.setName("out_critical_damping_pos_ori_port");
    out_critical_damping_pos_ori_port.doc("Output port for critical damping");
    out_critical_damping_pos_ori_port.setDataSample(dCriticalPosOri);
    ports()->addPort(out_critical_damping_pos_ori_port);

    out_endEffectorVelError_port.setName("out_endEffectorVelError_port");
    out_endEffectorVelError_port.doc("End Effector Velocity Error");
    out_endEffectorVelError_port.setDataSample(endEffectorVelError);
    ports()->addPort(out_endEffectorVelError_port);

    out_systemEnergy_port.setName("out_systemEnergy_port");
    out_systemEnergy_port.doc("Total System Energy");
    out_systemEnergy_port.setDataSample(total_system_energy);
    ports()->addPort(out_systemEnergy_port);

    out_desPos_curPos_port.setName("out_desPos_curPos_port");
    out_desPos_curPos_port.doc("DesPos and CurPos");
    out_desPos_curPos_port.setDataSample(desPos_curPos);
    ports()->addPort(out_desPos_curPos_port);

	out_Fctrl_port.setName("out_Fctrl_port");
    out_Fctrl_port.doc("Control Force");
    out_Fctrl_port.setDataSample(Fctrl_XYZ);
    ports()->addPort(out_Fctrl_port);


    
}
void CartImpedanceControl::displayCurrentState() {
    cc->displayCurrentState();
    dc->displayCurrentState();
    jd->displayCurrentState();
    return;
}

void CartImpedanceControl::setGains_anisotropic(const Eigen::VectorXf& Stiffness_Des, const Eigen::VectorXf& Damping_Des) {

    this->Stiffness_Des[0] = Stiffness_Des[0];
    this->Stiffness_Des[1] = Stiffness_Des[1];
    this->Stiffness_Des[2] = Stiffness_Des[2];
    this->Stiffness_Des[3] = Stiffness_Des[3];
    this->Stiffness_Des[4] = Stiffness_Des[4];
    this->Stiffness_Des[5] = Stiffness_Des[5];

    this->Damping_Des[0] = Damping_Des[0];
    this->Damping_Des[1] = Damping_Des[1];
    this->Damping_Des[2] = Damping_Des[2];
    this->Damping_Des[3] = Damping_Des[3];
    this->Damping_Des[4] = Damping_Des[4];
    this->Damping_Des[5] = Damping_Des[5];
}

void CartImpedanceControl::setGains(float gainP, float gainD, float orientation_gainP, float orientation_gainD) {
    assert(gainP >= 0);
    assert(gainD >= 0);
    assert(orientation_gainP >= 0);
    assert(orientation_gainD >= 0);
    this->gainP = gainP;
    this->gainD = gainD;
    this->orientation_gainP = orientation_gainP;
    this->orientation_gainD = orientation_gainD;



}

/// NEW PART OF THE CODE //
void CartImpedanceControl::quatToEuler(const double qW, const double qX, const double qY, const double qZ, Eigen::VectorXf &EulerAngle){

    double sinr_cosp, cosr_cosp, sinp, siny_cosp, cosy_cosp;
    double roll, pitch, yaw;

    sinr_cosp = 0; cosr_cosp = 0; sinp = 0; siny_cosp = 0; cosy_cosp = 0;
    roll = 0; pitch = 0; yaw = 0;

    // ROLL
    sinr_cosp = +2.0 * (qW * qX + qY * qZ);
    cosr_cosp = +1.0 - 2.0 * (qX * qX + qY * qY);
    roll = atan2(sinr_cosp, cosr_cosp);

    // PITCH
    sinp = +2.0 * (qW * qY - qZ * qX);
    if (fabs(sinp) >= 1)
        pitch = copysign(M_PI / 2, sinp); // use 90 degrees if out of range
    else
        pitch = asin(sinp);

    // YAW
    siny_cosp = +2.0 * (qW * qZ + qX * qY);
    cosy_cosp = +1.0 - 2.0 * (qY * qY + qZ * qZ);
    yaw = atan2(siny_cosp, cosy_cosp);

    EulerAngle[0] = roll;
    EulerAngle[1] = pitch;
    EulerAngle[2] = yaw;
}

void CartImpedanceControl::calculate_K_Max_Beta_Pos_Ori(Eigen::VectorXf &F_MAX, Eigen::VectorXf &Max_Disp, Eigen::VectorXf &Tau_MAX, Eigen::VectorXf &Max_Disp_Ori, Eigen::MatrixXf &K_Max_System, Eigen::VectorXf &Beta ,Eigen::MatrixXf &K_Max_System_Ori, Eigen::VectorXf &Beta_Ori, int i) {

    K_Max_System(i,i) = F_MAX[i]/Max_Disp[i];
	//K_Max_System(i,i) = 200;    
	Beta[i] = sqrt(log(K_Max_System(i,i))/pow(Max_Disp[i],2));
	if (std::isnan(Beta[i])){
		Beta[i] = 0;
	}

    K_Max_System_Ori(i,i) = Tau_MAX[i]/Max_Disp_Ori[i];
    Beta_Ori[i] = sqrt(log(K_Max_System_Ori(i,i))/pow(Max_Disp_Ori[i],2));
	if (std::isnan(Beta_Ori[i])){
		Beta_Ori[i] = 0;
	}

}

void CartImpedanceControl::get_K_Total_Pos_Ori(Eigen::VectorXf &Beta, Eigen::MatrixXf &K_Var_Pos, Eigen::MatrixXf &K_Total_Des_Pos, Eigen::VectorXf &Beta_Ori, Eigen::MatrixXf &K_Var_Ori, Eigen::MatrixXf &K_Total_Des_Ori, int i) {

    K_Var_Pos(i,i) = exp(pow(Beta[i] * (position_error[i]), 2));
    K_Total_Des_Pos(i,i) = kConstPos(i,i) + K_Var_Pos(i,i);

    K_Var_Ori(i, i) = exp(pow(Beta_Ori[i] * (angular_error[i]),2));
    K_Total_Des_Ori(i,i) = kConstOri(i,i) + K_Var_Ori(i,i);

}

void CartImpedanceControl::calculate_F_Max(Eigen::VectorXf &Max_Disp, Eigen::VectorXf &F_Max, int i){
		
		//F_Max[i] = 20;
		
		// CALIBRATED //
		if (i == 0){
		if (Max_Disp[i] >= 0.003){
			F_Max[i] = 30;		
		}else if(Max_Disp[i] >= 0.001 && Max_Disp[i] < 0.003) {
			F_Max[i] = 15;
		}else {
			F_Max[i] = 12.5;
		}
		}else if (i == 1){

		if (Max_Disp[i] >= 0.09){
			F_Max[i] = 30;		
		}else if(Max_Disp[i] >= 0.05 && Max_Disp[i] < 0.09) {
			F_Max[i] = 25;
		}else if(Max_Disp[i] >= 0.04 && Max_Disp[i] < 0.05) {
			F_Max[i] = 22.5;
		}else if(Max_Disp[i] >= 0.03 && Max_Disp[i] < 0.04) {
			F_Max[i] = 20;
		}else if(Max_Disp[i] >= 0.02 && Max_Disp[i] < 0.03) {
			F_Max[i] = 18;
		}else if(Max_Disp[i] >= 0.007 && Max_Disp[i] < 0.02) {
			F_Max[i] = 12.5;
		}else if(Max_Disp[i] >= 0.006 && Max_Disp[i] < 0.007) {
			F_Max[i] = 10;
		}else if(Max_Disp[i] >= 0.004 && Max_Disp[i] < 0.006) {
			F_Max[i] = 6;
		}else if(Max_Disp[i] >= 0.003 && Max_Disp[i] < 0.004) {
			F_Max[i] = 5;
		}else if(Max_Disp[i] >= 0.002 && Max_Disp[i] < 0.003) {
			F_Max[i] = 4;
		}else if(Max_Disp[i] >= 0.001 && Max_Disp[i] < 0.002) {
			F_Max[i] = 3;
		}else {
			F_Max[i] = 2;
		}
		}else{

		if (Max_Disp[i] >= 0.08){
			F_Max[i] = 30;		
		}else if(Max_Disp[i] >= 0.05 && Max_Disp[i] < 0.08) {
			F_Max[i] = 28; 
		}else if(Max_Disp[i] >= 0.04 && Max_Disp[i] < 0.05) {
			F_Max[i] = 26;
		}else if(Max_Disp[i] >= 0.03 && Max_Disp[i] < 0.04) {
			F_Max[i] = 24;
		}else if(Max_Disp[i] >= 0.02 && Max_Disp[i] < 0.03) {
			F_Max[i] = 20;
		}else if(Max_Disp[i] >= 0.007 && Max_Disp[i] < 0.02) {
			F_Max[i] = 12.5;
		}else if(Max_Disp[i] >= 0.006 && Max_Disp[i] < 0.007) {
			F_Max[i] = 10;
		}else if(Max_Disp[i] >= 0.004 && Max_Disp[i] < 0.006) {
			F_Max[i] = 6;
		}else if(Max_Disp[i] >= 0.003 && Max_Disp[i] < 0.004) {
			F_Max[i] = 5;
		}else if(Max_Disp[i] >= 0.002 && Max_Disp[i] < 0.003) {
			F_Max[i] = 4;
		}else if(Max_Disp[i] >= 0.001 && Max_Disp[i] < 0.002) {
			F_Max[i] = 3;
		}else {
			F_Max[i] = 2;
		}
}


	// NOT CALIBRATED //
	/*if (i == 0){
		if (Max_Disp[i] >= 0.003){
			F_Max[i] = 20;		
		}else if(Max_Disp[i] >= 0.002 && Max_Disp[i] < 0.003) {
			F_Max[i] = 18;
		}else if (Max_Disp[i] >= 0.001 && Max_Disp[i] < 0.002){
			F_Max[i] = 12;
		}else {
			F_Max[i] = 12;
		}
	
	}else{
		if (Max_Disp[i] >= 0.08){
        	slopeFmax = 50;
        	F_Max[i] = (-slopeFmax * (0.08-Max_Disp[i]) + 11)/1;

    	}else if(Max_Disp[i] >= 0.06 && Max_Disp[i] < 0.08){
        slopeFmax = 100;
        F_Max[i] = (-slopeFmax * (0.06-Max_Disp[i]) + 9)/1;

    }else if(Max_Disp[i] >= 0.04 && Max_Disp[i] < 0.06){
        slopeFmax = 50;
        F_Max[i] = (-slopeFmax * (0.04-Max_Disp[i]) + 8)/1;

    }else if(Max_Disp[i] >= 0.02 && Max_Disp[i] < 0.04){
        slopeFmax = 100;
        F_Max[i] = (-slopeFmax * (0.02-Max_Disp[i]) + 6)/1;

    }else if(Max_Disp[i] >= 0.01 && Max_Disp[i] < 0.02){
        slopeFmax = 100;
        F_Max[i] = (-slopeFmax * (0.01-Max_Disp[i]) + 5)/1;
		//F_Max[i] = 12;
		
    }else if(Max_Disp[i] >= 0.007 && Max_Disp[i] < 0.01){
        slopeFmax = 333;
        F_Max[i] = (-slopeFmax * (0.007-Max_Disp[i]) + 4)/1;
		
    }else if(Max_Disp[i] >= 0.006 && Max_Disp[i] < 0.007){
        slopeFmax = 1000;
        F_Max[i] = (-slopeFmax * (0.006-Max_Disp[i]) + 3)/1;

    }else if(Max_Disp[i] >= 0.001 && Max_Disp[i] < 0.006){
        slopeFmax = 250;
        F_Max[i] = -slopeFmax * (0.001-Max_Disp[i]) + 3.5;

    }else{
	    F_Max[i] = 1;
	}
}*/


}

void CartImpedanceControl::calculate_Tau_Max(Eigen::VectorXf &Max_Disp_Ori, Eigen::VectorXf &Tau_Max, int i){

	//Tau_Max[i] = 1;

	// CALIBRATED //
	if (i == 0){	
		Tau_Max[i] = 20;
		}
	else{
		if (Max_Disp_Ori[i] >= 0.261799){
			Tau_Max[i] = 20;
		}else if (Max_Disp_Ori[i] >= 0.0872665 && Max_Disp_Ori[i] < 0.261799){
			Tau_Max[i] = 15;
		}else {
			Tau_Max[i] = 5;
	}
  }



	// NOT CALIBRATED //
    /*if (Max_Disp_Ori[i] >= 0.08){
        slopeTauMax_Ori = 50;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.08-Max_Disp_Ori[i]) + 11)/2;

    }else if(Max_Disp_Ori[i] >= 0.06 && Max_Disp_Ori[i] < 0.08){
        slopeTauMax_Ori = 100;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.06-Max_Disp_Ori[i]) + 9)/2;

    }else if(Max_Disp_Ori[i] >= 0.04 && Max_Disp_Ori[i] < 0.06){
        slopeTauMax_Ori = 50;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.04-Max_Disp_Ori[i]) + 8)/2;

    }else if(Max_Disp_Ori[i] >= 0.02 && Max_Disp_Ori[i] < 0.04){
        slopeTauMax_Ori = 100;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.02-Max_Disp_Ori[i]) + 6)/2;

    }else if(Max_Disp_Ori[i] >= 0.01 && Max_Disp_Ori[i] < 0.02){
        slopeTauMax_Ori = 100;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.01-maxDisp_Ori[i]) + 5)/2;

    }else if(maxDisp_Ori[i] >= 0.007 && maxDisp_Ori[i] < 0.01){
        slopeTauMax_Ori = 333;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.007-maxDisp_Ori[i]) + 4)/2;

    }else if(maxDisp_Ori[i] >= 0.006 && maxDisp_Ori[i] < 0.007){
        slopeTauMax_Ori = 1000;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.006-maxDisp_Ori[i]) + 3)/2;

    }else if(maxDisp_Ori[i] >= 0.001 && maxDisp_Ori[i] < 0.006){
        slopeTauMax_Ori = 250;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.001-maxDisp_Ori[i]) + 1)/2;

    }else
        Tau_Max[i] = 1.5/2;*/

    /*if (Max_Disp_Ori[i] >= 0.08){
        slopeTauMax_Ori = 50;
        Tau_Max[i] =(-slopeTauMax_Ori * (0.08-Max_Disp_Ori[i]) + 11)/(1);

    }else if(Max_Disp_Ori[i] >= 0.06 && Max_Disp_Ori[i] < 0.08){
        slopeTauMax_Ori = 100;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.06-Max_Disp_Ori[i]) + 9)/(1);

    }else if(Max_Disp_Ori[i] >= 0.04 && Max_Disp_Ori[i] < 0.06){
        slopeTauMax_Ori = 50;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.04-Max_Disp_Ori[i]) + 8)/(1);

    }else if(Max_Disp_Ori[i] >= 0.02 && Max_Disp_Ori[i] < 0.04){
        slopeTauMax_Ori = 100;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.02-Max_Disp_Ori[i]) + 6)/(1);

    }else if(Max_Disp_Ori[i] >= 0.01 && Max_Disp_Ori[i] < 0.02){
        slopeTauMax_Ori = 100;
        Tau_Max[i] = 0.8; //0.0125(-slopeFmax_Ori * (0.01-maxDisp_Ori[i]) + 5)/30;

    }else if(maxDisp_Ori[i] >= 0.007 && maxDisp_Ori[i] < 0.01){
    slopeFmax_Ori = 333;
    Fmax_Ori[i] = (-slopeFmax_Ori * (0.007-maxDisp_Ori[i]) + 4)/2;

    }else if(maxDisp_Ori[i] >= 0.006 && maxDisp_Ori[i] < 0.007){
    slopeFmax_Ori = 1000;
    Fmax_Ori[i] = (-slopeFmax_Ori * (0.006-maxDisp_Ori[i]) + 3)/2;

    }else if(maxDisp_Ori[i] >= 0.001 && maxDisp_Ori[i] < 0.006){
    slopeFmax_Ori = 250;
    Fmax_Ori[i] = -slopeFmax_Ori * (0.001-maxDisp_Ori[i]) + 1;

    }else
        Tau_Max[i] = 0.0125;*/
}

void CartImpedanceControl::boxShaping_Online(Eigen::VectorXf &Max_Disp, int i){
	
	
		if (currentTime >= time0 + 1){
		//RTT::log(RTT::Error) << i << RTT::endlog() << RTT::endlog();	
		Max_Disp[i] = Max_Disp[i] - 0.0025;
		if (i == 2){
			time0 = currentTime;
		}    	
    	if(Max_Disp[i] <= 0.01) {
    		Max_Disp[i] = 0.01;
    	}
	}
}

void CartImpedanceControl::trajectory_execution(Eigen::VectorXf &DesPos, Eigen::VectorXf &EndEff_VelLinear, int i){

	if (cc->des_poses_var[0] != numDesPos_x || cc->des_poses_var[1] != numDesPos_y ||
        cc->des_poses_var[2] != numDesPos_z) {
        numDesPos_x = cc->des_poses_var[0];
        numDesPos_y = cc->des_poses_var[1];
        numDesPos_z = cc->des_poses_var[2];
    }

	if (numTraj_CircTraj_Act == 0) {
		if (i == 0){
            DesPos[0] = numDesPos_x;
            EndEff_VelLinear[0] = -cc->twists_var[0];
        } else if(i == 1){
            DesPos[1] = numDesPos_y + lineDispTrajCoef * sin(omegaTraj * 2 * 3.14 * (currentTime - trajStartTime));
            EndEff_VelLinear[1] = -cc->twists_var[1] + omegaTraj * 2 * 3.14 * lineDispTrajCoef * cos(omegaTraj * 2 * 3.14 * (currentTime - trajStartTime));
        }else{
            DesPos[2] = numDesPos_z;
            EndEff_VelLinear[2] = -cc->twists_var[2];
        }
	}else {
		if (i == 0) {
            DesPos[0] = numDesPos_x  + circRadiusTraj * cos(omegaTraj * 2 * 3.14 * (currentTime - trajStartTime));
            EndEff_VelLinear[0] = -cc->twists_var[0]  - omegaTraj * 2 * 3.14 * circRadiusTraj * sin(omegaTraj * 2 * 3.14 * (currentTime - trajStartTime));
		}else if(i == 1){
            DesPos[1] = numDesPos_y + circRadiusTraj * sin(omegaTraj * 2 * 3.14 * (currentTime - trajStartTime)) ;
            EndEff_VelLinear[1] = -cc->twists_var[1] + omegaTraj * 2 * 3.14 * circRadiusTraj * cos(omegaTraj * 2 * 3.14 * (currentTime - trajStartTime));
        } else {
            DesPos[2] = numDesPos_z ;
            EndEff_VelLinear[2] = -cc->twists_var[2];
        }
		}  
}

void CartImpedanceControl::calculation_K1_K2_Kd_Max_PMax_PMid_Max_Dist_Vec(Eigen::MatrixXf &K1, Eigen::MatrixXf &K2, Eigen::MatrixXf &Kd_Max, Eigen::VectorXf &Pmax, Eigen::VectorXf &Pmid, Eigen::VectorXf &Max_Dist_Vec, Eigen::VectorXf &Beta, Eigen::MatrixXf &K_Max_System, Eigen::VectorXf &Max_Disp, int i) {

    //if ((velLinearError[i] <= 0.005 && fabs(position_error[i]) > fabs(preEndEff_PosError3x1[i])) || fabs(position_error[i]) > fabs(preEndEff_PosError3x1[i])) {
if (fabs(position_error[i]) > fabs(preEndEff_PosError3x1[i])) {
        Max_Dist_Vec[i] = curPos3x1[i] - desPos3x1[i];
        Kd_Max(i,i) = kConstPos(i,i) + kVarPos(i, i);
    } else{

        Max_Dist_Vec[i] = Max_Dist_Vec[i];
        Kd_Max(i,i) = Kd_Max(i,i);
    }

    // DEFINING Pmax3x1 and Pmid3x1
    //if (curPos3x1[i] - desPos3x1[i] >= 0){

        Pmax[i] = desPos3x1[i] + Max_Dist_Vec[i];
        Pmid[i] = desPos3x1[i] + (Pmax[i] - desPos3x1[i])/2;

    /*}else{

        Pmax[i] = desPos3x1[i] + Max_Dist_Vec[i];
        Pmid[i] = desPos3x1[i] + (Pmax[i] - desPos3x1[i])/2;

    }*/

    // DEFINING K1_3x3 AND K2_3x3 w.r.t MAXIMUM STIFFNESS AT THE RELEASED POINT
	K1(i, i) = (2 / pow(Max_Dist_Vec[i], 2)) * ((exp(pow(Beta[i], 2) * pow(Max_Dist_Vec[i], 2)) - pow(Beta[i], 2) * pow(Max_Dist_Vec[i], 2) + pow(Beta[i], 2) * kConstPos(i, i) * pow(Max_Dist_Vec[i], 2)))/pow(Beta[i], 2);
    K2(i, i) = (2 / pow(Max_Dist_Vec[i], 2)) * ((exp(pow(Beta[i], 2) * pow(Max_Dist_Vec[i], 2)) - pow(Beta[i], 2) * pow(Max_Dist_Vec[i], 2) + pow(Beta[i], 2) * kConstPos(i, i) * pow(Max_Dist_Vec[i], 2)))/pow(Beta[i], 2);
    
    if (K1(i,i) >= kConstPos(i,i) + kVarPos(i, i)){
        K1(i,i) = kConstPos(i,i) + kVarPos(i, i);
    }

    if (K2(i,i) >= kConstPos(i,i) + kVarPos(i, i)){
        K2(i,i) = kConstPos(i,i) + kVarPos(i, i);
    }

    //K1(i, i) = 2 * Kd_Max(i, i);
    //K2(i, i) = 2 * Kd_Max(i, i);

}

void CartImpedanceControl::calculation_K1_K2_Kd_Max_PMax_PMid_Max_Dist_Vec_Ori(Eigen::MatrixXf &K1_Ori, Eigen::MatrixXf &K2_Ori, Eigen::MatrixXf &Kd_Max_Ori, Eigen::VectorXf &Pmax_Ori, Eigen::VectorXf &Pmid_Ori, Eigen::VectorXf &Max_Dist_Vec_Ori, Eigen::VectorXf &Beta_Ori, Eigen::MatrixXf &K_Max_System_Ori, int i) {

    //if ((velAngularError[i] <= 0.005 && fabs(angular_error[i]) > fabs(preEndEff_OriError3x1[i])) || fabs(angular_error[i]) > fabs(preEndEff_OriError3x1[i])) {
    if ( fabs(angular_error[i]) > fabs(preEndEff_OriError3x1[i])) {
        Max_Dist_Vec_Ori[i] = curOri3x1[i] - desOri3x1[i];
        Kd_Max_Ori(i,i) = kConstOri(i,i) + kVarOri(i, i);
    } else{

        Max_Dist_Vec_Ori[i] = Max_Dist_Vec_Ori[i];
        Kd_Max_Ori(i,i) = Kd_Max_Ori(i,i);
    }

    // DEFINING Pmax3x1 and Pmid3x1
    //if (curOri3x1[i] - desOri3x1[i] >= 0){

        Pmax_Ori[i] = desOri3x1[i] + Max_Dist_Vec_Ori[i];
        Pmid_Ori[i] = desOri3x1[i] + (Pmax_Ori[i] - desOri3x1[i])/2;
        //RTT::log(RTT::Error) << "Condition 01" << RTT::endlog() << RTT::endlog();
    /*}else{

        Pmax_Ori[i] = desOri3x1[i] + Max_Dist_Vec_Ori[i];
        Pmid_Ori[i] = desOri3x1[i] + (Pmax_Ori[i] - desOri3x1[i])/2;
        //RTT::log(RTT::Error) << "Condition 02" << RTT::endlog() << RTT::endlog();
    }*/

    // DEFINING K1_Ori_3x3 AND K2_Ori_3x3 w.r.t MAXIMUM STIFFNESS AT THE RELEASED POINT
	K1_Ori(i, i) = (2 / pow(Max_Dist_Vec_Ori[i], 2)) * ((exp(pow(Beta_Ori[i], 2) * pow(Max_Dist_Vec_Ori[i], 2)) - pow(Beta_Ori[i], 2) * pow(Max_Dist_Vec_Ori[i], 2) + pow(Beta_Ori[i], 2) * kConstOri(i, i) * pow(Max_Dist_Vec_Ori[i], 2)))/pow(Beta_Ori[i], 2);
    K2_Ori(i, i) = (2 / pow(Max_Dist_Vec_Ori[i], 2)) * ((exp(pow(Beta_Ori[i], 2) * pow(Max_Dist_Vec_Ori[i], 2)) - pow(Beta_Ori[i], 2) * pow(Max_Dist_Vec_Ori[i], 2) + pow(Beta_Ori[i], 2) * kConstOri(i, i) * pow(Max_Dist_Vec_Ori[i], 2)))/pow(Beta_Ori[i], 2);
   
	if (K1_Ori(i,i) >= kConstOri(i,i) + kVarOri(i, i)){
        K1_Ori(i,i) = kConstOri(i,i) + kVarOri(i, i);
    }

    if (K2_Ori(i,i) >= kConstOri(i,i) + kVarOri(i, i)){
        K2_Ori(i,i) = kConstOri(i,i) + kVarOri(i, i);
    } 
    //K1_Ori(i, i) = 2 * Kd_Max_Ori(i, i);
    //K2_Ori(i, i) = 2 * Kd_Max_Ori(i, i);
}

void CartImpedanceControl::calculation_Kc_DampConstLiner_FMidCtrl_FkCtrl(Eigen::MatrixXf &Kc, Eigen::MatrixXf &DampConstLinear, Eigen::VectorXf &FMidCtrl, Eigen::VectorXf &FkCtrl, int i) {


        
            //if((std::abs(position_error[i]) >= 0 && std::abs(velLinearError[i]) > 0)){// && preEndEff_VelError3x1[i] >= 0){
            if(sgn(position_error[i],0.005f) == sgn(velLinearError[i],0.005f)){
              Kc.setZero();
                DampConstLinear(i,i) = 0;
                FMidCtrl[i] = 0;

            }
//else if(){// && preEndEff_VelError3x1[i] <= 0){
            //    Kc_3x3.setZero();
            //    DampConstLinear(i,i) = 0;
            //    FMidCtrl[i] = 0;

            //}else{
            //    if((position_error[i] >= 0 /*&& preEndEff_VelError3x1[i] >= 0*/ && velLinearError[i] < 0 && curPos3x1[i] >= Pmid3x1[i]) || (position_error[i] <= 0 /*&& preEndEff_VelError3x1[i] <= 0*/ && velLinearError[i] > 0 && curPos3x1[i] <= Pmid3x1[i])){
            //        //RTT::log(RTT::Error) << "Condition 1" << RTT::endlog() << RTT::endlog();
            //        /*if (fabs(velLinearError[i]) < 0.2){
            //            DampConstLinear(i,i) = 0;
            //        }*/
            //        Kc = K2_3x3;
            //        FMidCtrl = Kc * (curPos3x1 - Pmid3x1);
            //        FkCtrl[i] = 0;

            //    }else if((position_error[i] >= 0 /*&& preEndEff_VelError3x1[i] >= 0*/ && velLinearError[i] < 0 && curPos3x1[i] < Pmid3x1[i]) || (position_error[i] <= 0 /*&& preEndEff_VelError3x1[i] <= 0*/ && velLinearError[i] > 0 && curPos3x1[i] > Pmid3x1[i])){
                    //RTT::log(RTT::Error) << "Condition 2" << RTT::endlog() << RTT::endlog();
                    /*if (fabs(velLinearError[i]) < 0.2){
                        DampConstLinear(i,i) = 0;
                    }*/
            //        Kc = K1_3x3;
            //        FMidCtrl = Kc * (Pmid3x1 - curPos3x1);
            //        FkCtrl[i] = 0;

            //    }
            //}
            else {//if((position_error[i] >= 0 && velLinearError[i] < 0) || (position_error[i] <= 0 && velLinearError[i] > 0)){
                    Kc = K2_3x3;
                    FMidCtrl = Kc * (Pmid3x1 - curPos3x1);
                    FkCtrl[i] = 0;
            }
        
    
}

void CartImpedanceControl::calculation_Kc_DampConstAngular_FMidCtrl_FkCtrl_Ori(Eigen::MatrixXf &Kc_Ori, Eigen::MatrixXf &DampConstAngular, Eigen::VectorXf &FMidCtrl_Ori, Eigen::VectorXf &FkCtrl_Ori, int i) {

    
            //if((std::abs(angular_error[i]) >= 0 && std::abs(velAngularError[i]) > 0)){// && preEndEff_VelError3x1[i] >= 0){
            if(sgn(angular_error[i]) == sgn(velAngularError[i])){
                Kc_Ori.setZero();
                DampConstAngular(i,i) = 0;
                FMidCtrl_Ori[i] = 0;

            }
//else if(){// && preEndEff_VelError3x1[i] <= 0){
            //    Kc_3x3.setZero();
            //    DampConstLinear(i,i) = 0;
            //    FMidCtrl[i] = 0;

            //}else{
            //    if((position_error[i] >= 0 /*&& preEndEff_VelError3x1[i] >= 0*/ && velLinearError[i] < 0 && curPos3x1[i] >= Pmid3x1[i]) || (position_error[i] <= 0 /*&& preEndEff_VelError3x1[i] <= 0*/ && velLinearError[i] > 0 && curPos3x1[i] <= Pmid3x1[i])){
            //        //RTT::log(RTT::Error) << "Condition 1" << RTT::endlog() << RTT::endlog();
            //        /*if (fabs(velLinearError[i]) < 0.2){
            //            DampConstLinear(i,i) = 0;
            //        }*/
            //        Kc = K2_3x3;
            //        FMidCtrl = Kc * (curPos3x1 - Pmid3x1);
            //        FkCtrl[i] = 0;

            //    }else if((position_error[i] >= 0 /*&& preEndEff_VelError3x1[i] >= 0*/ && velLinearError[i] < 0 && curPos3x1[i] < Pmid3x1[i]) || (position_error[i] <= 0 /*&& preEndEff_VelError3x1[i] <= 0*/ && velLinearError[i] > 0 && curPos3x1[i] > Pmid3x1[i])){
                    //RTT::log(RTT::Error) << "Condition 2" << RTT::endlog() << RTT::endlog();
                    /*if (fabs(velLinearError[i]) < 0.2){
                        DampConstLinear(i,i) = 0;
                    }*/
            //        Kc = K1_3x3;
            //        FMidCtrl = Kc * (Pmid3x1 - curPos3x1);
            //        FkCtrl[i] = 0;

            //    }
            //}
            else {//if((position_error[i] >= 0 && velLinearError[i] < 0) || (position_error[i] <= 0 && velLinearError[i] > 0)){
                    Kc_Ori = K2_3x3;
                    FMidCtrl_Ori = Kc_Ori * (Pmid3x1_Ori - curOri3x1);
                    FkCtrl_Ori[i] = 0;
            }
       
            /*if((angular_error[i]) >= 0 && velAngularError[i] >= 0 && preEndEff_VelOriError3x1[i] >= 0){
                Kc_Ori.setZero();
                DampConstAngular(i,i) = 0;
                FMidCtrl_Ori[i] = 0;

            }else if((angular_error[i]) <= 0 && velAngularError[i] <= 0 && preEndEff_VelOriError3x1[i] <= 0){
                Kc_Ori.setZero();
                DampConstAngular(i,i) = 0;
                FMidCtrl_Ori[i] = 0;

            }else{
                if((angular_error[i] >= 0 && preEndEff_VelOriError3x1[i] >= 0 && velAngularError[i] < 0 && curOri3x1[i] >= Pmid3x1_Ori[i]) || (angular_error[i] <= 0 && preEndEff_VelOriError3x1[i] <= 0 && velAngularError[i] > 0 && curOri3x1[i] <= Pmid3x1_Ori[i])){
                    //RTT::log(RTT::Error) << "Condition 3" << RTT::endlog() << RTT::endlog();
                    Kc_Ori = K2_Ori_3x3;
                    FMidCtrl_Ori = Kc_Ori * (curOri3x1 - Pmid3x1_Ori);
                    //F_midPointErrorCtrl_Ori[0] = F_midPointErrorCtrl_Ori[0]/100;
                    FkCtrl_Ori[i] = 0;

                }else if((angular_error[i] >= 0 && preEndEff_VelOriError3x1[i] >= 0 && velAngularError[i] < 0 && curOri3x1[i] < Pmid3x1_Ori[i]) || (angular_error[i] <= 0 && preEndEff_VelOriError3x1[i] <= 0 && velAngularError[i] > 0 && curOri3x1[i] > Pmid3x1_Ori[i])){
                    //RTT::log(RTT::Error) << "Condition 4" << RTT::endlog() << RTT::endlog();
                    Kc_Ori = K1_Ori_3x3;
                    FMidCtrl_Ori = Kc_Ori * (Pmid3x1_Ori - curOri3x1);
                    //F_midPointErrorCtrl_Ori[0] = F_midPointErrorCtrl_Ori[0]/100;
                    FkCtrl_Ori[i] = 0;

                }
            }*/
        
    
}

// // // // // // // // //

double CartImpedanceControl::getSimulationTime() {
    return 1E-9 * RTT::os::TimeService::ticks2nsecs(RTT::os::TimeService::Instance()->getTicks());

}

void CartImpedanceControl::compute() {

    Eigen::VectorXf error(6*cc->number_of_end_effectors);
    Eigen::VectorXf vel_error(6*cc->number_of_end_effectors);
    
    Eigen::MatrixXf inertia_c_inv = dc->in_inertia_var.inverse();
	lambda_c = pinv->compute(cc->in_jacobian_var * inertia_c_inv * cc->in_jacobian_var.transpose(), thresh);

    kConstPos(0,0) = Stiffness_Des[0]; kConstPos(1,1) = Stiffness_Des[1]; kConstPos(2,2) = Stiffness_Des[2];
    kConstOri(0,0) = Stiffness_Des[3]; kConstOri(1,1) = Stiffness_Des[4]; kConstOri(2,2) = Stiffness_Des[5];
    dConstLinear(0,0)= Damping_Des[0]; dConstLinear(1,1)= Damping_Des[1]; dConstLinear(2,2)= Damping_Des[2];
    dConstAngular(0,0)= Damping_Des[3]; dConstAngular(1,1)= Damping_Des[4]; dConstAngular(2,2)= Damping_Des[5];

	// NULLSPACE CTRL ELEMENTS - GAINS AND DESIRED JOINT VALUES
    float gainNullP, gainNullD;
    gainNullP = 100;
    gainNullD = 20;

    desJntPos7x1[0] = -0.154303;
    desJntPos7x1[1] = -0.300414;
    desJntPos7x1[2] = 0.222139;
    desJntPos7x1[3] = 1.21592;
    desJntPos7x1[4] = 0.0785831;
    desJntPos7x1[5] = 1.27652;
    desJntPos7x1[6] = 1.12973;

	// DELTA TIME CALCULATION + CURRENT TIME
    deltaTime = time2 - time1;
    currentTime=getSimulationTime() -timeInitial;

    for (int i = 0; i < cc->number_of_end_effectors; i++) {

        int pose_index = i * 7;

        if(cc->poses_var.segment<4>(pose_index+3).transpose()*cc->des_poses_var.segment<4>(pose_index+3)>0){

            cc->des_poses_var.segment<4>(pose_index+3) = -cc->des_poses_var.segment<4>(pose_index+3);
        }

		Eigen::Matrix3f skew;
        CosimaUtilities::skewMatrix(cc->poses_var.segment<3>(pose_index + 4), skew);
        angular_error = cc->des_poses_var(pose_index + 3) * cc->poses_var.segment<3>(pose_index + 4) - cc->poses_var(pose_index + 3) * cc->des_poses_var.segment<3>(pose_index + 4) - skew * cc->des_poses_var.segment<3>(pose_index + 4);

        vel_error.segment<6>(i*6) = cc->des_twists_var.segment<6>(i * 6) - cc->twists_var.segment<6>(i * 6);
        
		for (int i = 0; i < 3; i++) {

            // TRAJECTORY ACTIVATION
            if (trajectory_activation == 1){
                trajectory_execution(desPos3x1, velLinearError, i);
            }
		    //
		    curPos3x1[i] = cc->poses_var[i];

			// Desired and current end-effector Linear Velocity
			desVel3x1[i] = cc->des_twists_var[i];
			curVel3x1[i] = cc->twists_var[i];

		    // Desired and Current end-effector Orientation
		    //quatToEuler(cc->des_poses_var[3],cc->des_poses_var[4],cc->des_poses_var[5],cc->des_poses_var[6], desOri3x1);
		    //quatToEuler(cc->poses_var[3],cc->poses_var[4],cc->poses_var[5],cc->poses_var[6], curOri3x1);

        	position_error[i] = desPos3x1[i] - curPos3x1[i];
        	//angular_error[i] = desOri3x1[i] - curOri3x1[i];


		}

		
		//  RTT::log(RTT::Error) << "desPos3x1" << desPos3x1.transpose() << RTT::endlog() ;
		
        // Position and Orientation Error
        //position_error = cc->des_poses_var.segment<3>(pose_index) - cc->poses_var.segment<3>(pose_index);
        //angular_error = desOri3x1 - curOri3x1;	
		
		for (int i = 0; i < 3; i++){
        // Computation of end-effector linear/angular velocity errors
        velLinearError[i] = cc->des_twists_var[i] - cc->twists_var[i];
        velAngularError[i] = cc->des_twists_var[i + 3] - cc->twists_var[i + 3];		
		endEffectorVelError[i] = velLinearError[i];
        endEffectorVelError[i+3] = velAngularError[i];

		// Matrices of Constant Stiffness gains for end-effector Position and Orientation
        //kConstPos(i,i) = gainP;
        //kConstOri(i,i) = orientation_gainP;

        // Matrices of Constant Damping gains for end-effector Linear and Angular Velocities
        //dConstLinear(i,i) = gainD;
        //dConstAngular(i, i) = orientation_gainD;

        // ONLINE BOX SHAPING
        if(onlineBoxShapingSwitch == 0){

        //Fmax[i] = 15;
        //TauMax[i] = 0.0125;

        // FMAX AND TAU_MAX CALCULATION
        calculate_F_Max (maxDisp, Fmax, i);
        calculate_Tau_Max (maxDisp_Ori, TauMax, i);
		//RTT::log(RTT::Error) << "DEBUGGING!" << RTT::endlog() << RTT::endlog();

        } else
        {
        	boxShaping_Online (maxDisp, i);
			calculate_F_Max (maxDisp, Fmax, i);
            calculate_Tau_Max (maxDisp_Ori, TauMax, i);

        }

        // K_SYSTEM_MAX, MAXIMUM DISPLACEMENT AND BETA CALCULATION FOR BOTH END-EEFECTOR POSITION AND ORIENTATION
        calculate_K_Max_Beta_Pos_Ori(Fmax, maxDisp, TauMax, maxDisp_Ori, KmaxSystem, beta, KmaxSystem_Ori, beta_Ori, i);

        // K_TOTAL AND K_VAR CALCULATION FOR BOTH POSITION AND OREINTATION
        get_K_Total_Pos_Ori(beta, kVarPos, kTotalDesPos, beta_Ori, kVarOri, kTotalDesOri, i);

		if (kTotalDesPos(i,i) > KmaxSystem(i,i) && beta[i]>0){
        	kVarPos(i,i) = Fmax[i]/fabs(position_error[i]) - kConstPos(i,i);
		}else if(kTotalDesPos(i,i) > KmaxSystem(i,i)){
                kVarPos(i,i) = 0;
        }

        if (kTotalDesOri(i,i) > KmaxSystem_Ori(i,i)&& beta_Ori[i]>0){
           kVarOri(i,i) = TauMax[i]/fabs(angular_error[i]) - kConstOri(i,i);
        }else if(kTotalDesOri(i,i) > KmaxSystem_Ori(i,i)){
                kVarOri(i,i) = 0;
        }

        //kTotal6x6(i,i) = kConstPos(i,i) + kVarPos(i,i);
        //kTotal6x6(i+3,i+3) = kConstOri(i,i) + kVarOri(i,i);
        dTotal6x6(i,i) = dConstLinear(i,i) +  dVarPos(i,i);
        dTotal6x6(i+3,i+3) = dConstAngular(i,i) + dVarOri(i,i);

		// TESTING - CRITICAL DAMPING CALCULATION TESTING
		//dCritical6x6 = 2 * (lambda_c * kTotal6x6).sqrt();		

	    // ATTRACTIVE FORCE FOR POSITION AND ORIENTATION CTRL CALCULATION
        Fk[i] = (kConstPos(i,i) + kVarPos(i,i)) * (position_error[i]);
        Fk_Ori[i] = (kConstOri(i,i) + kVarOri(i, i)) * angular_error[i];

        // DAMPING FORCES FOR LINEAR AND ANGULAR VELOCITIES - DEFINITIONS
        Fd_Pos[i] = dConstLinear(i,i) * velLinearError[i];
        Fd_Ori[i] = dConstAngular(i,i) * velAngularError[i];

        // TESTING - DAMPING FORCES FOR LINEAR AND ANGULAR VELOCITIES - USE OF CRITICAL DAMPING
        /*Fd_Pos[i] = dCritical6x6(i,i) * velLinearError[i];
        Fd_Ori[i] = dCritical6x6(i+3,i+3) * velAngularError[i];*/

		//////////////////////// PROPOSED METHOD - END EFFECTOR POSITION ////////////////////////
        calculation_K1_K2_Kd_Max_PMax_PMid_Max_Dist_Vec(K1_3x3, K2_3x3, KdMax, Pmax3x1, Pmid3x1, maxDistVec_3x1, beta, KmaxSystem, maxDisp, i);

        //////////////////////// PROPOSED METHOD - END EFFECTOR ORIENTATION ////////////////////////
        calculation_K1_K2_Kd_Max_PMax_PMid_Max_Dist_Vec_Ori(K1_Ori_3x3, K2_Ori_3x3, KdMax_Ori, Pmax3x1_Ori, Pmid3x1_Ori, maxDistVec_Ori_3x1, beta_Ori, KmaxSystem_Ori, i);

        ////////////////////////////////////////////////////////////////////////////////////////////
		
		//////////////////////// DEFINING THE CONDITIONS FOR THE PROPOSED METHOD - END EFFECTOR POSITION ////////////////////////
        calculation_Kc_DampConstLiner_FMidCtrl_FkCtrl(Kc_3x3, dConstLinear, F_midPointErrorCtrl, Fk, i);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //////////////////////// DEFINING THE CONDITIONS FOR THE PROPOSED METHOD - END EFFECTOR ORIENTATION ////////////////////////
        calculation_Kc_DampConstAngular_FMidCtrl_FkCtrl_Ori(Kc_Ori_3x3, dConstAngular, F_midPointErrorCtrl_Ori, Fk_Ori, i);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		


            // NEW PART OF THE CODE // 30.01.2019
			/*if (fabs(curVel3x1[i]) > fabs(desVel3x1[i])){
				if(position_error[i] > 0){
					Fdamping[i] = - 1 * dConstLinear(i,i) * fabs(velLinearError[i]);				
				}else if(position_error[i] < 0){
				if(position_error[i] > 0){
					Fdamping[i] =  1 * dConstLinear(i,i) * fabs(velLinearError[i]);				
				}
			}else{
				Fdamping[i] = 0;
				}
			}*/
            /////////////////////////////////////
////////////////////////////////////////////////////



////////////////////////////////////////////////////
////////////////
            }

        // NULLSPACE DESIRED JOINT VALUES AND CTRL GAINS INITIALIZATION
        for (int i = 0; i < 7; i++) {
            // Current Nullspace Joint Position and Velocity
            curJntPos7x1[i] = jd->in_robotstatus_var.angles[i];
            curJntVel7x1[i] = jd->in_robotstatus_var.velocities[i];

            // Matrices of Constant Stiffness and Damping gains for Nullspace Ctrl
            kConstNull (i,i) = gainNullP;
            dConstNull (i,i) = gainNullD;

        }

	nullSpaceCtrl = nullSpace->calculateNullspace(cc->in_jacobian_var) * ((kConstNull + kVarNull) * (desJntPos7x1 - curJntPos7x1) + (dConstNull + dVarNull) * (desJntVel7x1 - curJntVel7x1));


//////////////////////// NEW C++ LINES ADDED IN THIS SECTION ////////////////////////

        //////////////////////// DEFINING THE CONDITIONS FOR THE PROPOSED METHOD - END EFFECTOR POSITION ////////////////////////
           //if(currentTime>= 7){
               
        
        error.segment<3>((i * 6)) = Fk + F_midPointErrorCtrl;
        error.segment<3>((i*6)+3) = Fk_Ori + F_midPointErrorCtrl_Ori;

        vel_error.segment<3>((i*6)) = Fd_Pos;
		vel_error.segment<3>((i*6)+3) = Fd_Ori;

		// NULLSPACE CONTROL VARIABLE STIFFNESS
        //nullSpaceCtrl = nullSpace->calculateNullspace(cc->in_jacobian_var) * (kConstNull * (desJntPos7x1 - curJntPos7x1) /*- psi * x_t*/ + dConstNull * (desJntVel7x1 - curJntVel7x1));
		
	/*RTT::log(RTT::Error) << "Pmax3x1 = " << Pmax3x1.transpose() << RTT::endlog() << "Pmid3x1 = " << Pmid3x1.transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "KmaxSystem ==> " << KmaxSystem.diagonal().transpose() << RTT::endlog() << "beta ==> " << beta.transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "KdMax = " << (KdMax.diagonal()).transpose() << RTT::endlog();
    RTT::log(RTT::Error) << "K1_3x3 = " << (K1_3x3.diagonal()).transpose() << RTT::endlog() << "K2_3x3 = " << (K2_3x3.diagonal()).transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "Fmax ==> " << Fmax.transpose() << RTT::endlog() << "maxDisp ==> " << maxDisp.transpose() << RTT::endlog();
    RTT::log(RTT::Error) << "dConstLinear_3x3 ==> " << (dConstLinear.diagonal()).transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "Fk ==> " << Fk.transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "F_midPointErrorCtrl ==> " << F_midPointErrorCtrl.transpose() << RTT::endlog() << RTT::endlog();*/

    /*RTT::log(RTT::Error) << "Pmax3x1_Ori = " << Pmax3x1_Ori.transpose() << RTT::endlog() << "Pmid3x1_Ori = " << Pmid3x1_Ori.transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "KmaxSystem_Ori ==> " << KmaxSystem_Ori.diagonal().transpose() << RTT::endlog() << "beta_Ori ==> " << beta_Ori.transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "KdMax_Ori = " << (KdMax_Ori.diagonal()).transpose() << RTT::endlog();
    RTT::log(RTT::Error) << "K1_Ori_3x3 = " << (K1_Ori_3x3.diagonal()).transpose() << RTT::endlog() << "K2_Ori_3x3 = " << (K2_Ori_3x3.diagonal()).transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "TauMax ==> " << TauMax.transpose() << RTT::endlog() << "maxDisp_Ori ==> " << maxDisp_Ori.transpose() << RTT::endlog();
    RTT::log(RTT::Error) << "dConstAngular_3x3 ==> " << (dConstAngular.diagonal()).transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "Fk_Ori ==> " << Fk_Ori.transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "F_midPointErrorCtrl_Ori ==> " << F_midPointErrorCtrl_Ori.transpose() << RTT::endlog() << RTT::endlog();
		
	RTT::log(RTT::Error) << "Total System Energy ==> " << total_system_energy[0] << RTT::endlog() << RTT::endlog() ;*/

	//RTT::log(RTT::Error) << "maxDisp ==> " << maxDisp.transpose() << RTT::endlog() << RTT::endlog(); 

	
		
		
    }

	// System Total Energy

        for (int i = 0; i < 3; i++) {


            if (sgn(position_error[i], 0.005f) == sgn(velLinearError[i], 0.005f)) {
                posError_oriError_6x1[i] = position_error[i];
            }else {
                posError_oriError_6x1[i] = curPos3x1[i] - Pmid3x1[i];
            }

            if(sgn(angular_error[i]) == sgn(velAngularError[i])){
                posError_oriError_6x1[i+3] = angular_error[i];
            }else{
                posError_oriError_6x1[i+3] = curOri3x1[i] - Pmid3x1_Ori[i];
            }
        }

    for(int i = 0; i < 3; i++) {
        if (sgn(position_error[i], 0.000f) == sgn(velLinearError[i], 0.000f)) {
            kTotal6x6(i,i) = kConstPos(i,i) + kVarPos(i,i);
        } else {
            kTotal6x6(i,i) = Kc_3x3(i,i);
        }
        if(sgn(angular_error[i]) == sgn(velAngularError[i])){
            kTotal6x6(i+3,i+3) = kConstOri(i,i) + kVarOri(i,i);
        }else{
            kTotal6x6(i+3,i+3) = Kc_Ori_3x3(i,i);
        }
    }

	potentialEnergy = 0.5 * posError_oriError_6x1.transpose() * kTotal6x6 * posError_oriError_6x1;
    kineticEnergy = 0.5 * endEffectorVelError.transpose() * lambda_c * endEffectorVelError;
    energy_vector[0] = potentialEnergy + kineticEnergy;

	// REQUIRED CTRL JOINT TORQUES CALCULATION AND IMPLEMENTATON
    Eigen::MatrixXf JMcinv = cc->in_jacobian_var * inertia_c_inv;
	Eigen::VectorXf F = ((error +  vel_error) + lambda_c * (dc->in_acceleration_var + JMcinv * (dc->in_coriolis_var) - dc->in_jacobian_dot_var * jd->in_robotstatus_var.velocities));

    cc->out_torques_var.torques = cc->in_jacobian_var.transpose() * F + nullSpaceCtrl;
	trqSat->saturateTorqueRate(cc->out_torques_var.torques, estimated_torque_var);
    cc->out_torques_var.torques = estimated_torque_var;
    out_estimated_torque_port.write(estimated_torque_var);

	//RTT::log(RTT::Error) << "estimated_torque_var ==> " << estimated_torque_var.transpose()  << RTT::endlog() << RTT::endlog();
	
	
	// // // // // // // // // // // // // // //
    
    //////////////////////// NEW C++ LINES ADDED IN THIS SECTION ////////////////////////
    preEndEff_PosError3x1 = position_error;
    preEndEff_VelError3x1 = velLinearError;

    preEndEff_OriError3x1 = angular_error;
    preEndEff_VelOriError3x1 = velAngularError;
    //////////////////////// //////////////////////// ////////////////////////

    out_endEffectorVelError_port.write(endEffectorVelError);

    // make a vector of variable stiffness Pos and Ori
    varStiffVec.head(3) = (kConstPos + kVarPos).diagonal() ;
    varStiffVec.tail(3) = (kConstOri + kVarOri).diagonal() ;
    out_variable_stiffness_pos_ori_port.write(varStiffVec);
    
    // make a vector of Pos Error and Ori Error
    //posError_oriError_6x1.head(3) = position_error;
    //posError_oriError_6x1.tail(3) = angular_error;
    out_posError_oriError_port.write(posError_oriError_6x1);
	
	
	total_system_energy.head(3) = energy_vector;
	out_systemEnergy_port.write(total_system_energy);

	desPos_curPos.head(3) = desPos3x1;
	desPos_curPos.tail(3) = curPos3x1;
	out_desPos_curPos_port.write(desPos_curPos);

	Fctrl_XYZ = error;
	out_Fctrl_port.write(Fctrl_XYZ);

}

bool CartImpedanceControl::configureHook() {
    return cc->configureHook() && dc->configureHook() && jd->configureHook();
}
bool CartImpedanceControl::startHook() {
	timeInitial = getSimulationTime();
    time1 = getSimulationTime()-timeInitial;

    return cc->startHook() && dc->startHook() && jd->startHook();
}
void CartImpedanceControl::updateHook() {
    
	time2 = getSimulationTime()-timeInitial;

	if (!cc->getData() || !dc->getData() || !jd->getData()) {
        RTT::log(RTT::Error) << this->getName() << " Cart:" << std::to_string(cc->getData()) << "|| DC:"
                             << std::to_string(dc->getData()) << "|| JD:" << std::to_string(jd->getData())
                             << RTT::endlog();
        return;
    }

	desPos3x1 = cc->des_poses_var.head<3>();
    compute();

    //cc->out_torques_var.torques = dc->in_gravity_var;
    cc->out_torques_port.write(cc->out_torques_var);

    time1 = time2;
}

void CartImpedanceControl::stopHook() {
    cc->stopHook();
    dc->stopHook();
    jd->stopHook();
}

void CartImpedanceControl::cleanupHook() {
    cc->cleanupHook();
    dc->cleanupHook();
    jd->cleanupHook();
}
ORO_CREATE_COMPONENT_LIBRARY()
ORO_LIST_COMPONENT_TYPE(CartImpedanceControl)

// make a vector for desTraj and curPos Error
//posError_desTraj_curPos_6x1.head(3) = desTraj3x1;
//posError_desTraj_curPos_6x1.tail(3) = curPos3x1;
//out_posError_desTraj_curPos_port.write(posError_desTraj_curPos_6x1);

 // PRINT DATA FOR MATLAB ANALYSIS TEST (FOR FUTURE TXT FORMAT)
        /*RTT::log(RTT::Error) << "," << desPos3x1[0] << "," << curPos3x1[0] <<
                                "," << desPos3x1[1] << "," << curPos3x1[1] <<
                                "," << desPos3x1[2] << "," << curPos3x1[2] <<
                                "," << cc->des_twists_var[0]  << "," << cc->twists_var[0] <<
                                "," << cc->des_twists_var[1]  << "," << cc->twists_var[1] <<
                                "," << cc->des_twists_var[2]  << "," << cc->twists_var[2] <<RTT::endlog();*/;
