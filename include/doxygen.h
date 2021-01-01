/*! \file doxygen.h
 *  \brief Main page for Doygen documentation
 */
 
/** \mainpage Matrix-Solver Community Project
 *
 *
 * The purpose of this project is to construct a pool of problems and solvers for 
 * these problems. 
 * It is a community of matrices and solvers; hereafter we refer it as the Community. 
 * For a problem (the Community contains problem descriptions), there are a bunch of 
 * solvers (the Community contains algorithm descriptions and source codes) for this 
 * problem. A solver could be applied to different problems (the Community compares 
 * performance and convergence of the solver for various of problems)
 *
 * 
 * \section sec_multigrid Our goal and strategy 
 *
 * We can think of the Community in a Multigrid + Capitalism point of view.
 *
 * \subsection ss_fine Fine grid (free market stage)
 *
 * (1) Collect problems and solvers. Allow duplications, for example same solution 
 * algorithm, but different implementation or different programming languages. 
 *
 * (2) Do not enforce much regulation. Allow the market to be FREE. 
 * 
 * (3) Keep all the record: problem description, solver code, test results, etc. 
 * Do not throw anything away.
 *
 * (4) We are currently at the initial stage. We try to find a minimal set of 
 * standard or rules. And then we let the market to evolve freely.
 *
 *
 * \subsection ss_coarse Coarse grid (state capitalism stage)
 *
 * (1)	As the market evolves, we might find at certain time when the market is 
 * out-of-control. This basically means the "free market" is very successful. And 
 * now we need to give more restrict standard or regulation.
 *
 * (2) We still keep everything before standardized; just a new branch containing 
 * standardized solvers. Different versions marked with different colors for 
 * convenience.
 * 
 * (3) At certain stage, we might write professional level codes for the solvers 
 * and form a package. We are now working toward this goal. 
 *
 * (4)	Host an annual user-developer meeting and set up a solver competition. Start 
 * a journal for computational softwares. 
 *
 *
 * 
 * \section sec_guide General guidelines 
 *
 * 1. For the moment, this project is still restricted to our own group only. 
 * It does not have to be user-friendly to people outside of the group. But at 
 * least, for people in the group, the problems should be easy to understand 
 * and the solvers should be easy to compile and run. 
 *
 * 2. Write as much comments as possible for easy maintenance and usage. 
 * See \ref page_comment for more detailed instructions on how to add comment 
 * using Doxygen. 
 *
 * 3. We mainly interested in the problems with algebraic information only 
 * (maybe with some mesh information and/or local stiffness matrices). See 
 * \ref page_steps for solving steps.  
 * 
 * 4. We are mathematicians, not software engineers. We are going to do it as
 * mathematicians without worrying about the implementation quality too much. 
 * We don't need to build an industry quality software package. It is for 
 * academic (research and education) purpose. Version control is still very 
 * important for our collaborative work; see \ref page_comment for CVS.
 *
 * 5. All source codes are protected. If a user is interested in a particular 
 * solver, he/she needs to contact the author (could some one in our own group). 
 *
 */
 
 
/**
 * \page page_steps Steps to solve problems
 *
 * To solve a problem (some times we only have the coefficent matrix and the 
 * right-hand side), we can do the following five steps: 
 *
 * \section step1 Step 1: obtain a problem
 *
 *	 We can assemble directly if we have the mesh and PDE information; otherwise,
 * we read from disk files. 
 *  
 * \section step2 Step 2: check matrix property
 *
 * We run a sequence of tests to see: whether this matrix is symmetric, positive 
 * definite, sparse, etc. 
 *   
 * \section step3 Step 3: select a solver
 *
 * Using some artificial intelligence, we pick a solver which is suitable for 
 * the problem we are dealing with. This is the brain of this whole project. 
 * Prof. Xu will work on this part. 
 *   
 * \section step4 Step 4: solve the system
 *
 * Once the solver has been chosen, we solve the system. 
 *   
 * \section step5 Step 5: post processing
 *   
 * Graphical or text output.
 *
 */
 
 
/**
 * \page page_algorithms Algorithms
 *
 * Here we write all the algorithms we implemented and tested. 
 *
 */
 
 
/**
 * \page page_performance Performance
 *
 * Here we collect the numerical examples.
 *
 *
 
 
/**
 * \page page_comment Comment and Version Control
 *
 * \section sec_comment Adding comment
 *
 * We use Doxygen as our automatically documentation generator which will make our 
 * future maintainance minimized. You can obtain the software (Windows, Linux and 
 * OS X) as well as its manual on the official website 
 * 
 * http://www.doxygen.org
 *
 * For an oridinary user, Doxygen is completely trivial to use. We only need to use 
 * some special marker in the usual comment as we put in c-files. 
 * See test.c and matvecio.c for examples.  
 *
 * \section sec_cvs Version control
 *
 * In order to coordinate our collaboration, we use CVS as our version control tool. 
 * If you are willing to contribute to the project, you can send 
 * Chensong Zhang multigrid@me.com an email.
 *
 * CVS is easy to use and it is available at all platforms. For Linux and OS X users, 
 * set up the following environment variables.
 *
 * For bash/sh users:
 *
 * $ export CVSROOT = solver@nfs1.screms.math.psu.edu:/homelab/solver/cvsroot
 * 
 * $ export CVS_RSH=`type -p ssh`
 *
 * For csh/tcsh
 *
 * $ setenv CVSROOT solver@nfs1.screms.math.psu.edu:/homelab/solver/cvsroot
 * 
 * $ setenv CVS_RSH `which ssh`
 *
 * Windows version is supposed to be fool-proof and we don't include an instruction 
 * for Windows version here. 
 *
 * We use SSH for security. But we don't use passwd access, instead everybody should 
 * send me you public SSH authorization code. If you already have one, it should locate 
 * in your home directory at
 *
 * ~/.ssh/id_rsa.pub
 *
 * If you don't have one in .ssh, you can create one by 
 * 
 * $ ssh-keygen -t rsa
 *
 * The program will ask you a few questions, just hit on RETURN for every question. 
 *
 * If you don't even have .ssh direction, you need to create on 
 * 
 * $ mkdir ~/.ssh
 * 
 * When you find this ***.pub file, please send it to me. So I can give you access to 
 * the CVS server without passwd. Now we are all set to go. 
 *
 * Now just do 
 * 
 * $ cvs co msc 
 * 
 * will give you the current version of the package. 
 * 
 */