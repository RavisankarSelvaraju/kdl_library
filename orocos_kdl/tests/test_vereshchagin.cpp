// #include "solvertest.hpp"
#include "test_vereshchagin.hpp"
#include <frames_io.hpp>
#include <framevel_io.hpp>
#include <kinfam_io.hpp>
#include <random>
#include <time.h>
#include <utilities/utility.h>
#include <iostream>

using namespace KDL;


// void setUp()
// {
//     srand( (unsigned)time( NULL ));

//     chain1.addSegment(Segment("Segment 1", Joint("Joint 1", Joint::RotZ),
//                               Frame(Vector(0.0,0.0,0.0))));
//     chain1.addSegment(Segment("Segment 2", Joint("Joint 2", Joint::RotX),
//                               Frame(Vector(0.0,0.0,0.9))));
//     chain1.addSegment(Segment("Segment 3", Joint("Joint 3", Joint::None),
//                               Frame(Vector(-0.4,0.0,0.0))));
//     chain1.addSegment(Segment("Segment 4", Joint("Joint 4", Joint::RotX),
//                               Frame(Vector(0.0,0.0,1.2))));
//     chain1.addSegment(Segment("Segment 5", Joint("Joint 5", Joint::None),
//                               Frame(Vector(0.4,0.0,0.0))));
//     chain1.addSegment(Segment("Segment 6", Joint("Joint 6", Joint::RotZ),
//                               Frame(Vector(0.0,0.0,1.4))));
//     chain1.addSegment(Segment("Segment 7", Joint("Joint 7", Joint::RotX),
//                               Frame(Vector(0.0,0.0,0.0))));
//     chain1.addSegment(Segment("Segment 8", Joint("Joint 8", Joint::RotZ),
//                               Frame(Vector(0.0,0.0,0.4))));
//     chain1.addSegment(Segment("Segment 9", Joint("Joint 9", Joint::None),
//                               Frame(Vector(0.0,0.0,0.0))));

//     chain2.addSegment(Segment("Segment 1", Joint("Joint 1", Joint::RotZ),
//                               Frame(Vector(0.0,0.0,0.5))));
//     chain2.addSegment(Segment("Segment 2", Joint("Joint 2", Joint::RotX),
//                               Frame(Vector(0.0,0.0,0.4))));
//     chain2.addSegment(Segment("Segment 3", Joint("Joint 3", Joint::RotX),
//                               Frame(Vector(0.0,0.0,0.3))));
//     chain2.addSegment(Segment("Segment 4", Joint("Joint 4", Joint::RotX),
//                               Frame(Vector(0.0,0.0,0.2))));
//     chain2.addSegment(Segment("Segment 5", Joint("Joint 5", Joint::RotZ),
//                               Frame(Vector(0.0,0.0,0.1))));


//     chain3.addSegment(Segment("Segment 1", Joint("Joint 1", Joint::RotZ),
//                               Frame(Vector(0.0,0.0,0.0))));
//     chain3.addSegment(Segment("Segment 2", Joint("Joint 2", Joint::RotX),
//                               Frame(Vector(0.0,0.0,0.9))));
//     chain3.addSegment(Segment("Segment 3", Joint("Joint 3", Joint::RotZ),
//                               Frame(Vector(-0.4,0.0,0.0))));
//     chain3.addSegment(Segment("Segment 4", Joint("Joint 4", Joint::RotX),
//                               Frame(Vector(0.0,0.0,1.2))));
//     chain3.addSegment(Segment("Segment 5", Joint("Joint 5", Joint::None),
//                               Frame(Vector(0.4,0.0,0.0))));
//     chain3.addSegment(Segment("Segment 6", Joint("Joint 6", Joint::RotZ),
//                               Frame(Vector(0.0,0.0,1.4))));
//     chain3.addSegment(Segment("Segment 7", Joint("Joint 7", Joint::RotX),
//                               Frame(Vector(0.0,0.0,0.0))));
//     chain3.addSegment(Segment("Segment 8", Joint("Joint 8", Joint::RotZ),
//                               Frame(Vector(0.0,0.0,0.4))));
//     chain3.addSegment(Segment("Segment 9", Joint("Joint 9", Joint::RotY),
//                               Frame(Vector(0.0,0.0,0.0))));


//     chain4.addSegment(Segment("Segment 1", Joint("Joint 1", Vector(10,0,0), Vector(1,0,1),Joint::RotAxis),
//                               Frame(Vector(0.0,0.0,0.5))));
//     chain4.addSegment(Segment("Segment 2", Joint("Joint 2", Vector(0,5,0), Vector(1,0,0),Joint::RotAxis),
//                               Frame(Vector(0.0,0.0,0.4))));
//     chain4.addSegment(Segment("Segment 3", Joint("Joint 3", Vector(0,0,5), Vector(1,0,4),Joint::RotAxis),
//                               Frame(Vector(0.0,0.0,0.3))));
//     chain4.addSegment(Segment("Segment 4", Joint("Joint 4", Vector(0,0,0), Vector(1,0,0),Joint::RotAxis),
//                               Frame(Vector(0.0,0.0,0.2))));
//     chain4.addSegment(Segment("Segment 5", Joint("Joint 5", Vector(0,0,0), Vector(0,0,1),Joint::RotAxis),
//                               Frame(Vector(0.0,0.0,0.1))));



//     //chain definition for vereshchagin's dynamic solver
//     Joint rotJoint0 = Joint(Joint::RotZ, 1, 0, 0.01);
//     Joint rotJoint1 = Joint(Joint::RotZ, 1, 0, 0.01);

//     Frame refFrame(Rotation::RPY(0.0, 0.0, 0.0), Vector(0.0, 0.0, 0.0));
//     Frame frame1(Rotation::RPY(0.0, 0.0, 0.0), Vector(0.0, -0.4, 0.0));
//     Frame frame2(Rotation::RPY(0.0, 0.0, 0.0), Vector(0.0, -0.4, 0.0));

//     //chain segments
//     Segment segment1 = Segment(rotJoint0, frame1);
//     Segment segment2 = Segment(rotJoint1, frame2);

//     //rotational inertia around symmetry axis of rotation
//     RotationalInertia rotInerSeg1(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

//     //spatial inertia
//     RigidBodyInertia inerSegment1(0.3, Vector(0.0, -0.4, 0.0), rotInerSeg1);
//     RigidBodyInertia inerSegment2(0.3, Vector(0.0, -0.4, 0.0), rotInerSeg1);
//     segment1.setInertia(inerSegment1);
//     segment2.setInertia(inerSegment2);

//     //chain
//     chaindyn.addSegment(segment1);
//     chaindyn.addSegment(segment2);

// 	// Motoman SIA10 Chain (for IK singular value tests)
// 	motomansia10.addSegment(Segment(Joint(Joint::None),
// 									Frame::DH_Craig1989(0.0, 0.0, 0.36, 0.0)));
// 	motomansia10.addSegment(Segment(Joint(Joint::RotZ),
// 									Frame::DH_Craig1989(0.0, PI_2, 0.0, 0.0)));
// 	motomansia10.addSegment(Segment(Joint(Joint::RotZ),
// 									Frame::DH_Craig1989(0.0, -PI_2, 0.36, 0.0)));
// 	motomansia10.addSegment(Segment(Joint(Joint::RotZ),
// 									Frame::DH_Craig1989(0.0, PI_2, 0.0, 0.0)));
// 	motomansia10.addSegment(Segment(Joint(Joint::RotZ),
// 									Frame::DH_Craig1989(0.0, -PI_2, 0.36, 0.0)));
// 	motomansia10.addSegment(Segment(Joint(Joint::RotZ),
// 									Frame::DH_Craig1989(0.0, PI_2, 0.0, 0.0)));
// 	motomansia10.addSegment(Segment(Joint(Joint::RotZ),
// 									Frame::DH_Craig1989(0.0, -PI_2, 0.0, 0.0)));
// 	motomansia10.addSegment(Segment(Joint(Joint::RotZ),
// 									Frame(Rotation::Identity(),Vector(0.0,0.0,0.155))));

//     // Motoman SIA10 Chain with Mass Parameters (for forward dynamics tests)

//     //  effective motor inertia is included as joint inertia
//     static const double scale=1;
//     static const double offset=0;
//     static const double inertiamotorA=5.0;      // effective motor inertia kg-m^2
//     static const double inertiamotorB=3.0;      // effective motor inertia kg-m^2
//     static const double inertiamotorC=1.0;      // effective motor inertia kg-m^2
//     static const double damping=0;
//     static const double stiffness=0;

//     //  Segment Inertias
//     KDL::RigidBodyInertia inert1(15.0, KDL::Vector(0.0, -0.02, 0.0),                       // mass, CM
//                                  KDL::RotationalInertia(0.1, 0.05, 0.1, 0.0, 0.0, 0.0));   // inertia
//     KDL::RigidBodyInertia inert2(5.0, KDL::Vector(0.0, -0.02, -0.1),
//                                  KDL::RotationalInertia(0.01, 0.1, 0.1, 0.0, 0.0, 0.0));
//     KDL::RigidBodyInertia inert3(5.0, KDL::Vector(0.0, -0.05, 0.02),
//                                  KDL::RotationalInertia(0.05, 0.01, 0.05, 0.0, 0.0, 0.0));
//     KDL::RigidBodyInertia inert4(3.0, KDL::Vector(0.0, 0.02, -0.15),
//                                  KDL::RotationalInertia(0.1, 0.1, 0.01, 0.0, 0.0, 0.0));
//     KDL::RigidBodyInertia inert5(3.0, KDL::Vector(0.0, -0.05, 0.01),
//                                  KDL::RotationalInertia(0.02, 0.01, 0.02, 0.0, 0.0, 0.0));
//     KDL::RigidBodyInertia inert6(3.0, KDL::Vector(0.0, -0.01, -0.1),
//                                  KDL::RotationalInertia(0.1, 0.1, 0.01, 0.0, 0.0, 0.0));
//     KDL::RigidBodyInertia inert7(1.0, KDL::Vector(0.0, 0.0, 0.05),
//                                  KDL::RotationalInertia(0.01, 0.01, 0.1, 0.0, 0.0, 0.0));

//     motomansia10dyn.addSegment(Segment(Joint(Joint::RotZ, scale, offset, inertiamotorA, damping, stiffness),
//                                Frame::DH(0.0, PI_2, 0.36, 0.0),
//                                inert1));
//     motomansia10dyn.addSegment(Segment(Joint(Joint::RotZ, scale, offset, inertiamotorA, damping, stiffness),
//                                Frame::DH(0.0, -PI_2, 0.0, 0.0),
//                                inert2));
//     motomansia10dyn.addSegment(Segment(Joint(Joint::RotZ, scale, offset, inertiamotorB, damping, stiffness),
//                                Frame::DH(0.0, PI_2, 0.36, 0.0),
//                                inert3));
//     motomansia10dyn.addSegment(Segment(Joint(Joint::RotZ, scale, offset, inertiamotorB, damping, stiffness),
//                                Frame::DH(0.0, -PI_2, 0.0, 0.0),
//                                inert4));
//     motomansia10dyn.addSegment(Segment(Joint(Joint::RotZ, scale, offset, inertiamotorC, damping, stiffness),
//                                Frame::DH(0.0, PI_2, 0.36, 0.0),
//                                inert5));
//     motomansia10dyn.addSegment(Segment(Joint(Joint::RotZ, scale, offset, inertiamotorC, damping, stiffness),
//                                Frame::DH(0.0, -PI_2, 0.0, 0.0),
//                                inert6));
//     motomansia10dyn.addSegment(Segment(Joint(Joint::RotZ, scale, offset, inertiamotorC, damping, stiffness),
//                                Frame::DH(0.0, 0.0, 0.0, 0.0),
//                                inert7));
//     motomansia10dyn.addSegment(Segment(Joint(Joint::None),
//                                        Frame(Rotation::Identity(),Vector(0.0,0.0,0.155))));

//     /** 
//      * KUKA LWR 4 Chain with Dynamics Parameters (for Forward Dynamics and Vereshchagin solver tests)
//      * Necessary test model for the Vereshchagin solver: KDL's implementation of the Vereshchagin solver 
//      * can only work with the robot chains that have equal number of joints and segments.
//      * Note: Joint effective inertia values in this model are closely aligned with the joint inertia
//      * of the real robot. These parameters are published in: Jubien, A., Gautier, M. and Janot, A., 
//      * "Dynamic identification of the Kuka LWR robot using motor torques and joint torque sensors data.",
//      * IFAC Proceedings Volumes, 2014., 47(3), pp.8391-8396.
//      */
// 	kukaLWR.addSegment(Segment(Joint(Joint::RotZ, scale, offset, 3.19, damping, stiffness),
// 				  Frame::DH_Craig1989(0.0, 1.5707963, 0.0, 0.0),
// 				  Frame::DH_Craig1989(0.0, 1.5707963, 0.0, 0.0).Inverse()*RigidBodyInertia(2,
// 												 Vector::Zero(),
// 												 RotationalInertia(0.0,0.0,0.0115343,0.0,0.0,0.0))));

// 	kukaLWR.addSegment(Segment(Joint(Joint::RotZ, scale, offset, 3.05, damping, stiffness),
// 				  Frame::DH_Craig1989(0.0, -1.5707963, 0.4, 0.0),
// 				  Frame::DH_Craig1989(0.0, -1.5707963, 0.4, 0.0).Inverse()*RigidBodyInertia(2,
// 												   Vector(0.0,-0.3120511,-0.0038871),
// 												   RotationalInertia(-0.5471572,-0.0000302,-0.5423253,0.0,0.0,0.0018828))));

// 	kukaLWR.addSegment(Segment(Joint(Joint::RotZ, scale, offset, 1.98, damping, stiffness),
// 				  Frame::DH_Craig1989(0.0, -1.5707963, 0.0, 0.0),
// 				  Frame::DH_Craig1989(0.0, -1.5707963, 0.0, 0.0).Inverse()*RigidBodyInertia(2,
// 												   Vector(0.0,-0.0015515,0.0),
// 												   RotationalInertia(0.0063507,0.0,0.0107804,0.0,0.0,-0.0005147))));

// 	kukaLWR.addSegment(Segment(Joint(Joint::RotZ, scale, offset, 2.05, damping, stiffness),
// 				  Frame::DH_Craig1989(0.0, 1.5707963, 0.39, 0.0),
// 				  Frame::DH_Craig1989(0.0, 1.5707963, 0.39, 0.0).Inverse()*RigidBodyInertia(2,
// 												   Vector(0.0,0.5216809,0.0),
// 												   RotationalInertia(-1.0436952,0.0,-1.0392780,0.0,0.0,0.0005324))));

// 	kukaLWR.addSegment(Segment(Joint(Joint::RotZ, scale, offset, 0.787, damping, stiffness),
// 				  Frame::DH_Craig1989(0.0, 1.5707963, 0.0, 0.0),
// 				  Frame::DH_Craig1989(0.0, 1.5707963, 0.0, 0.0).Inverse()*RigidBodyInertia(2,
// 												   Vector(0.0,0.0119891,0.0),
// 												   RotationalInertia(0.0036654,0.0,0.0060429,0.0,0.0,0.0004226))));

// 	kukaLWR.addSegment(Segment(Joint(Joint::RotZ, scale, offset, 0.391, damping, stiffness),
// 				  Frame::DH_Craig1989(0.0, -1.5707963, 0.0, 0.0),
// 				  Frame::DH_Craig1989(0.0, -1.5707963, 0.0, 0.0).Inverse()*RigidBodyInertia(2,
// 												   Vector(0.0,0.0080787,0.0),
// 												   RotationalInertia(0.0010431,0.0,0.0036376,0.0,0.0,0.0000101))));

// 	kukaLWR.addSegment(Segment(Joint(Joint::RotZ, scale, offset, 0.394, damping, stiffness),
// 				               Frame::Identity(),
// 				               RigidBodyInertia(2, Vector::Zero(), RotationalInertia(0.000001,0.0,0.0001203,0.0,0.0,0.0))));
// }



// void VereshchaginTest()
// {
//     std::cout << "KDL Vereshchagin Hybrid Dynamics Tests" <<std::endl;

//     // ########################################################################################
//     // Vereshchagin solver test 1
//     // ########################################################################################
//     /**
//      * Compute Hybrid Dynamics for KUKA LWR 4.
//      * 
//      * Test setup:
//      * - Operational-space task imposes acceleration constraints on the end-effector
//      * - External forces and feedforward joint torques are acting on the robot's structure, 
//      *   as disturbances from the environment
//      * 
//      * Expected result:
//      * - The solver computes the required joint torque commands (joint constraint torques)
//      *   that satisfy imposed acceleration constraints and at the same time, compensate for 
//      *   the above mentioned disturbances
//      * 
//      * Method to evaluate:
//      * - Compare the _resultant_ Cartesian accelerations of the end-effector's segment with 
//      *   the task-specified acceleration constraints
//      */
//     int solver_return = 0;
//     double eps = 1.e-3;

//     unsigned int nj = kukaLWR.getNrOfJoints();
//     unsigned int ns = kukaLWR.getNrOfSegments();

//     // Necessary test for the used robot model: KDL's implementation of the Vereshchagin solver 
//     // can only work with the robot chains that have equal number of joints and segments
//     // CPPUNIT_ASSERT(Equal(nj, ns));
//     // CPPUNIT_ASSERT(Equal(nj, ns));
//     CPPUNIT_ASSERT(Equal(nj, ns));
//     // CPPUNIT_ASSERT(Equal(nj, ns));
//     // CPPUNIT_ASSERT(Equal(nj, ns));

//     // Joint position, velocity, acceleration, feed-forward and constraint torques
//     KDL::JntArray q(nj); //input
//     KDL::JntArray qd(nj); //input
//     KDL::JntArray qdd(nj); //output
//     KDL::JntArray ff_tau(nj); //input
//     KDL::JntArray constraint_tau(nj); //output

//     // Random configuration
//     q(0) =  1.6;
//     q(1) =  0.0;
//     q(2) = -1.6;
//     q(3) = -1.57;
//     q(4) =  0.0;
//     q(5) =  1.57;
//     q(6) = -0.8;

//     qd(0) =  1.0;
//     qd(1) = -2.0;
//     qd(2) =  3.0;
//     qd(3) = -4.0;
//     qd(4) =  5.0;
//     qd(5) = -6.0;
//     qd(6) =  7.0;

//     // Random feedforwad torques acting on robot joints
//     ff_tau(0) =  5.0;
//     ff_tau(1) =  0.0;
//     ff_tau(2) =  0.0;
//     ff_tau(3) =  0.0;
//     ff_tau(4) =  0.0;
//     ff_tau(5) = -6.0;
//     ff_tau(6) =  0.0;

//     // External Wrench acting on the end-effector, expressed in base link coordinates
//     // Vereshchagin solver expects that external wrenches are expressed w.r.t. robot's base frame
//     KDL::Vector f(10.0, 15.0, 0.0);
//     KDL::Vector n(0.0, 0.0, 5.0);
//     KDL::Wrench f_tool(f, n);
//     KDL::Wrenches f_ext(ns);
//     f_ext[ns - 1] = f_tool; //input

//     /**
//      * Definition of the Cartesian Acceleration Constraints imposed on the end-effector.
//      * Note: the Vereshchagin solver expects that the input values in alpha parameters 
//      * (unit constraint forces) are expressed w.r.t. robot's base frame. 
//      * However, the acceleration energy setpoints, i.e. beta parameters, are expressed w.r.t. above
//      * defined unit constraint forces. More specifically, each DOF (element) in beta parameter corresponds
//      * to its respective DOF (column) of the unit constraint force matrix (alpha).
//     */
//     int number_of_constraints = 6;

//     // Constraint Unit forces defined for the end-effector
//     Jacobian alpha_unit_force(number_of_constraints);

//     // Set directions in which the constraint force should work. Alpha in the solver
//     Twist unit_force_x_l(
//         Vector(1.0, 0.0, 0.0), 
//         Vector(0.0, 0.0, 0.0));
//     alpha_unit_force.setColumn(0, unit_force_x_l); // constraint active

//     Twist unit_force_y_l(
//         Vector(0.0, 1.0, 0.0),
//         Vector(0.0, 0.0, 0.0));
//     alpha_unit_force.setColumn(1, unit_force_y_l); // constraint active

//     Twist unit_force_z_l(
//         Vector(0.0, 0.0, 1.0),
//         Vector(0.0, 0.0, 0.0));
//     alpha_unit_force.setColumn(2, unit_force_z_l); // constraint active

//     Twist unit_force_x_a(
//         Vector(0.0, 0.0, 0.0),
//         Vector(0.0, 0.0, 0.0));
//     alpha_unit_force.setColumn(3, unit_force_x_a); // constraint disabled... In this direction, end-effector's motion is left to emerge naturally

//     Twist unit_force_y_a(
//         Vector(0.0, 0.0, 0.0),
//         Vector(0.0, 0.0, 0.0));
//     alpha_unit_force.setColumn(4, unit_force_y_a); // constraint disabled... In this direction, end-effector's motion is left to emerge naturally

//     Twist unit_force_z_a(
//         Vector(0.0, 0.0, 0.0),
//         Vector(0.0, 0.0, 1.0));
//     alpha_unit_force.setColumn(5, unit_force_z_a); // constraint active

//     // Acceleration energy for the end-effector.
//     JntArray beta_energy(number_of_constraints);
//     beta_energy(0) = -0.5;
//     beta_energy(1) = -0.5;
//     beta_energy(2) =  0.0;
//     beta_energy(3) =  0.0; // this value has no impact on computations, since its corresponding constraint is disabled
//     beta_energy(4) =  0.0; // this value has no impact on computations, since its corresponding constraint is disabled
//     beta_energy(5) =  0.2;

//     // Arm root acceleration (robot's base mounted on an even surface)
//     // Note: Vereshchagin solver takes root acc. with opposite sign comparead to the KDL's FD and RNE solvers
//     Twist root_Acc(Vector(0.0, 0.0, 9.81), Vector(0.0, 0.0, 0.0));

//     ChainHdSolver_Vereshchagin vereshchaginSolver(kukaLWR, root_Acc, number_of_constraints);
//     solver_return = vereshchaginSolver.CartToJnt(q, qd, qdd, alpha_unit_force, beta_energy, f_ext, ff_tau, constraint_tau);
//     if (solver_return < 0) std::cout << "KDL: Vereshchagin solver ERROR: " << solver_return << std::endl;

//     // ########################################################################################
//     // Final comparison of the _resultant_ end-effector's Cartesian accelerations
//     // and the task-specified acceleration constraints

//     // Number of frames on the robot = ns + 1
//     std::vector<Twist> xDotdot(ns + 1);
//     // This solver's function returns Cartesian accelerations of links in robot base coordinates
//     vereshchaginSolver.getTransformedLinkAcceleration(xDotdot);
  
//     // Additional getters for the intermediate solver's outputs: Useful for state- simulation and estimation purposes
//     // Magnitude of the constraint forces acting on the end-effector: Lagrange Multiplier
//     Eigen::VectorXd nu(number_of_constraints);
//     vereshchaginSolver.getContraintForceMagnitude(nu);
   
   
//     // Total torque acting on each joint
//     JntArray total_tau(nj);
//     vereshchaginSolver.getTotalTorque(total_tau);
 
  

//     // ########################################################################################
//     // Vereshchagin solver test 2
//     // ########################################################################################
//     Vector constrainXLinear(1.0, 0.0, 0.0);
//     Vector constrainXAngular(0.0, 0.0, 0.0);
//     Vector constrainYLinear(0.0, 0.0, 0.0);
//     Vector constrainYAngular(0.0, 0.0, 0.0);
//     // Vector constrainZLinear(0.0, 0.0, 0.0);
//     //Vector constrainZAngular(0.0, 0.0, 0.0);
//     Twist constraintForcesX(constrainXLinear, constrainXAngular);
//     Twist constraintForcesY(constrainYLinear, constrainYAngular);
//     //Twist constraintForcesZ(constrainZLinear, constrainZAngular);
//     Jacobian alpha(1);
//     //alpha.setColumn(0, constraintForcesX);
//     alpha.setColumn(0, constraintForcesX);
//     //alpha.setColumn(0, constraintForcesZ);

//     //Acceleration energy at  the end-effector
//     JntArray betha(1); //set to zero
//     betha(0) = 0.0;
//     //betha(1) = 0.0;
//     //betha(2) = 0.0;

//     //arm root acceleration
//     Vector linearAcc(0.0, 10, 0.0); //gravitational acceleration along Y
//     Vector angularAcc(0.0, 0.0, 0.0);
//     Twist twist1(linearAcc, angularAcc);

//     //external forces on the arm
//     Vector externalForce1(0.0, 0.0, 0.0);
//     Vector externalTorque1(0.0, 0.0, 0.0);
//     Vector externalForce2(0.0, 0.0, 0.0);
//     Vector externalTorque2(0.0, 0.0, 0.0);
//     Wrench externalNetForce1(externalForce1, externalTorque1);
//     Wrench externalNetForce2(externalForce2, externalTorque2);
//     Wrenches externalNetForce;
//     externalNetForce.push_back(externalNetForce1);
//     externalNetForce.push_back(externalNetForce2);
//     //~Definition of constraints and external disturbances
//     //-------------------------------------------------------------------------------------//


//     //Definition of solver and initial configuration
//     //-------------------------------------------------------------------------------------//
//     int numberOfConstraints = 1;
//     ChainHdSolver_Vereshchagin constraintSolver(chaindyn, twist1, numberOfConstraints);

//     //These arrays of joint values contain actual and desired values
//     //actual is generated by a solver and integrator
//     //desired is given by an interpolator
//     //error is the difference between desired-actual
//     //in this test only the actual values are printed.
//     const int k = 1;
//     JntArray jointPoses[k];
//     JntArray jointRates[k];
//     JntArray jointAccelerations[k];
//     JntArray jointFFTorques[k];
//     JntArray jointConstraintTorques[k];
//     for (int i = 0; i < k; i++)
//     {
//         JntArray jointValues(chaindyn.getNrOfJoints());
//         jointPoses[i] = jointValues;
//         jointRates[i] = jointValues;
//         jointAccelerations[i] = jointValues;
//         jointFFTorques[i] = jointValues;
//         jointConstraintTorques[i] = jointValues;
//     }

//     // Initial arm position configuration/constraint
//     JntArray jointInitialPose(chaindyn.getNrOfJoints());
//     jointInitialPose(0) = 0.0; // initial joint0 pose
//     jointInitialPose(1) = PI/6.0; //initial joint1 pose, negative in clockwise
//     //j0=0.0, j1=pi/6.0 correspond to x = 0.2, y = -0.7464
//     //j0=2*pi/3.0, j1=pi/4.0 correspond to x = 0.44992, y = 0.58636

//     //actual
//     jointPoses[0](0) = jointInitialPose(0);
//     jointPoses[0](1) = jointInitialPose(1);

//     //~Definition of solver and initial configuration
//     //-------------------------------------------------------------------------------------//


//     //Definition of process main loop
//     //-------------------------------------------------------------------------------------//
//     //Time required to complete the task move(frameinitialPose, framefinalPose)
//     double taskTimeConstant = 0.1;
//     double simulationTime = 1*taskTimeConstant;
//     double timeDelta = 0.01;

//     const std::string msg = "Assertion failed, check matrix and array sizes";

//     for (double t = 0.0; t <=simulationTime; t = t + timeDelta)
//     {
//         // CPPUNIT_ASSERT_EQUAL((int)SolverI::E_NOERROR, constraintSolver.CartToJnt(jointPoses[0], jointRates[0], jointAccelerations[0], alpha, betha, externalNetForce, jointFFTorques[0], jointConstraintTorques[0]));
//         // CPPUNIT_ASSERT_EQUAL((int)SolverI::E_NOERROR, constraintSolver.CartToJnt(jointPoses[0], jointRates[0], jointAccelerations[0], alpha, betha, externalNetForce, jointFFTorques[0], jointConstraintTorques[0]));
//         CPPUNIT_ASSERT_EQUAL((int)SolverI::E_NOERROR, constraintSolver.CartToJnt(jointPoses[0], jointRates[0], jointAccelerations[0], alpha, betha, externalNetForce, jointFFTorques[0], jointConstraintTorques[0]));
//         // CPPUNIT_ASSERT_EQUAL((int)SolverI::E_NOERROR, constraintSolver.CartToJnt(jointPoses[0], jointRates[0], jointAccelerations[0], alpha, betha, externalNetForce, jointFFTorques[0], jointConstraintTorques[0]));
//         // CPPUNIT_ASSERT_EQUAL((int)SolverI::E_NOERROR, constraintSolver.CartToJnt(jointPoses[0], jointRates[0], jointAccelerations[0], alpha, betha, externalNetForce, jointFFTorques[0], jointConstraintTorques[0]));

//         //Integration(robot joint values for rates and poses; actual) at the given "instanteneous" interval for joint position and velocity.
//         jointRates[0](0) = jointRates[0](0) + jointAccelerations[0](0) * timeDelta; //Euler Forward
//         jointPoses[0](0) = jointPoses[0](0) + (jointRates[0](0) - jointAccelerations[0](0) * timeDelta / 2.0) * timeDelta; //Trapezoidal rule
//         jointRates[0](1) = jointRates[0](1) + jointAccelerations[0](1) * timeDelta; //Euler Forward
//         jointPoses[0](1) = jointPoses[0](1) + (jointRates[0](1) - jointAccelerations[0](1) * timeDelta / 2.0) * timeDelta;
//         jointFFTorques[0] = jointConstraintTorques[0];
//         //printf("time, j0_pose, j1_pose, j0_rate, j1_rate, j0_acc, j1_acc, j0_constraintTau, j1_constraintTau \n");
//         printf("%f          %f      %f       %f     %f       %f      %f     %f      %f\n", t, jointPoses[0](0), jointPoses[0](1), jointRates[0](0), jointRates[0](1), jointAccelerations[0](0), jointAccelerations[0](1), jointConstraintTorques[0](0), jointConstraintTorques[0](1));
//     }
// }

using namespace std;
int main(int argc, char **argv)
{
    cout<<"Hello World"<<endl;
}