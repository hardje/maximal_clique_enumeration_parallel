MCE
Maximal Clique Enumeration

Current version: 
	Version 1.0
	August 23, 2016 
	
Complete Documentation can be found in mceDOCs.txt
	
	
Description:
	This project consists of three files.
	
	adjacency.cu:
		Takes in two text files of real number data and generates 
		an adjacency matrix file from them.
	mce.cu: 
		Takes in an adjacency matrix file and determines all of the
		maximal cliques that exist within the data.
	mcMeasure.cu:
		Takes in a set of maximal cliques and calculates a measure.
	
Libraries:	
	<iostream>: 	 Standard input output functions
	<fstream>:		 Read from and write to files	
	<sstream>:		 Provides functions for working with strings
	<string>: 		 Provides access to string data type
	<stdlib.h>: 	 Provides use of size_t data type
	<thrust/scan.h>: Provides parallel prefix scan function
	<vector>: 		 Provides access to the vector data structure
		
	
Set up and operation:
	Requirements:
		Compute Capability 3.0 or higher GPU device
    
    To compile:
		Requires use of Nvidia's CUDA Compiler: NVCC
			http://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/#axzz4FS98r9VS
        Go to file directory containing the cu files and type command:
        
		nvcc -arch=sm_35 -rdc=true <filename.cu> -o executableName
		
		to compile each cu file needed.
			
	Running:
		Once all of the files needed have been compiled, use the following parameters to run them.
		Command Line Parameters:
		
		adjacency.cu:
			MANDATORY PARAMETERS:
				-f [complete path to file1] [complete path to file2]
					Specify the input data files
				-o [int > 0]
					Specify the number of value per feature vector
			OPTIONAL PARAMETERS:
				-e [real number >= 0] 
					Specify the Epsilon value to be used when finding neighbours
					(Euclidean Distance measure is used, two vectors must be within Epsilon to be neighbours)
					Defaults to 0.4
				
		mce.cu:
			MANDATORY PARAMETERS:
				-v [int > 0]: The number of vertices in the input file
			OPTIONAL PARAMETERS:
				-gpu [int >= 0]: Specify which GPU device to run device code on
			 		(must be a valid device id)
			 		Defaults to device 0
 
 				-f [full path and file name]:
 					Specifies input file.
 					Defaults to adjMatrix.txt
 
 				-help: Prints out all command parameters, then exits program.
 				
 				-sc: Specify that cliques consisting of only one node should not
 						written to the output file
		
		mcMeasure.cu:
			MANDATORY PARAMETERS:
				-v [int > 0]: The number of vertices in the input file
			OPTIONAL PARAMETERS:
				-f [full path and file name]:
 					Specifies input file.
 					Defaults to cliques.txt
 				-sc: Specifies that cliques consisting of only one node
 					SHOULD factor into the measure calculation
 					(Note: This will cause the final measure to differ)
		
	Input Files:
		adjacency.cu:
			This program will take in exactly 2 files as input
			These files must adhere to the following rules:
				-Each input file is made up of real number values (Program will read them as floats)
				-Each input file is of the same dimensions as all other input files
				-Each file is laid out so that features of a feature vector are in order and adjecent
					to eachother in the file
					
					For example, for feature vectors made up of some features, X Y Z, the file should look like
						X0 Y0 Z0 X1 Y1 Z1 ... XN YN ZN  (Where N is the number of feature vectors)
					As opposed to:
						X0 X1 ... XN Y0 Y1 ... YN Z0 Z1 ... ZN
						
		mce.cu: 
			This program will take in an adjacency matrix text file as input.
			Such a file may be generated by adjacency.cu.
			This file must consist of only 1's and 0's (no spaces)
			The number of columns must be even and equal to the number of rows
			
		mcMeasure.cu:
			This file takes in a text file of maximal cliques as input.
			A maximal clique text file may be generated by mce.cu.
			Each row in the file consists of integer numbers between 
				0 and [(number of vectors in corresponding adjacency matrix) - 1]
				where these numbers represent the id's of vertices in the adjacency matrix.
			
		
Copyright and Licensing:
License	 : Licensed under the Non-Profit Open Software License version 3.0
 1) Grant of Copyright License. Licensor grants You a worldwide, royalty-free, 
 non-exclusive, sublicensable license, for the duration of the copyright, to do the following:

 a) to reproduce the Original Work in copies, either alone or as part of a collective work;

 b) to translate, adapt, alter, transform, modify, or arrange the Original Work, thereby 
 creating derivative works ("Derivative Works") based upon the Original Work;

 c) to distribute or communicate copies of the Original Work and Derivative Works 
 to the public, with the proviso that copies of Original Work or Derivative Works 
 that You distribute or communicate shall be licensed under this Non-Profit Open Software 
 License or as provided in section 17(d);

 d) to perform the Original Work publicly; and

 e) to display the Original Work publicly.

 2) Grant of Patent License. Licensor grants You a worldwide, royalty-free, 
 non-exclusive, sublicensable license, under patent claims owned or controlled by
  the Licensor that are embodied in the Original Work as furnished by the Licensor, 
  for the duration of the patents, to make, use, sell, offer for sale, have made,
   and import the Original Work and Derivative Works.

 3) Grant of Source Code License. The term "Source Code" means the preferred 
 form of the Original Work for making modifications to it and all available 
 documentation describing how to modify the Original Work. Licensor agrees to 
 provide a machine-readable copy of the Source Code of the Original Work along 
 with each copy of the Original Work that Licensor distributes. Licensor reserves 
 the right to satisfy this obligation by placing a machine-readable copy of the 
 Source Code in an information repository reasonably calculated to permit 
 inexpensive and convenient access by You for as long as Licensor continues 
 to distribute the Original Work.

 4) Exclusions From License Grant. Neither the names of Licensor, nor the names 
 of any contributors to the Original Work, nor any of their trademarks or service 
 marks, may be used to endorse or promote products derived from this Original Work 
 without express prior permission of the Licensor. Except as expressly stated
  herein, nothing in this License grants any license to Licensor's trademarks,
  copyrights, patents, trade secrets or any other intellectual property. No patent 
  license is granted to make, use, sell, offer for sale, have made, or import embodiments 
  of any patent claims other than the licensed claims defined in Section 2. No license 
 is granted to the trademarks of Licensor even if such marks are included in the Original
  Work. Nothing in this License shall be interpreted to prohibit Licensor from licensing 
  under terms different from this License any Original Work that Licensor otherwise would 
  have a right to license.

 5) External Deployment. The term "External Deployment" means the use, distribution, or
 communication of the Original Work or Derivative Works in any way such that the Original 
 Work or Derivative Works may be used by anyone other than You, whether those works are 
 distributed or communicated to those persons or made available as an application intended 
 for use over a network. As an express condition for the grants of license hereunder, 
 You must treat any External Deployment by You of the Original Work or a Derivative 
 Work as a distribution under section 1(c).

 6) Attribution Rights. You must retain, in the Source Code of any Derivative Works 
 that You create, all copyright, patent, or trademark notices from the Source Code of 
 the Original Work, as well as any notices of licensing and any descriptive text 
 identified therein as an "Attribution Notice." You must cause the Source Code for 
 any Derivative Works that You create to carry a prominent Attribution Notice reasonably 
 calculated to inform recipients that You have modified the Original Work.

 7) Warranty of Provenance and Disclaimer of Warranty. The Original Work is provided 
 under this License on an "AS IS" BASIS and WITHOUT WARRANTY, either express or implied, 
 including, without limitation, the warranties of non-infringement, merchantability or 
 fitness for a particular purpose. THE ENTIRE RISK AS TO THE QUALITY OF THE ORIGINAL WORK 
 IS WITH YOU. This DISCLAIMER OF WARRANTY constitutes an essential part of this License. 
 No license to the Original Work is granted by this License except under this disclaimer.

 8) Limitation of Liability. Under no circumstances and under no legal theory, whether 
 in tort (including negligence), contract, or otherwise, shall the Licensor be liable 
 to anyone for any direct, indirect, special, incidental, or consequential damages of 
 any character arising as a result of this License or the use of the Original Work 
 including, without limitation, damages for loss of goodwill, work stoppage, computer 
 failure or malfunction, or any and all other commercial damages or losses. This limitation 
 of liability shall not apply to the extent applicable law prohibits such limitation.

 9) Acceptance and Termination. If, at any time, You expressly assented to this License, 
 that assent indicates your clear and irrevocable acceptance of this License and all of 
 its terms and conditions. If You distribute or communicate copies of the Original Work 
 or a Derivative Work, You must make a reasonable effort under the circumstances to obtain 
 the express assent of recipients to the terms of this License. This License conditions 
 your rights to undertake the activities listed in Section 1, including your right to create 
 Derivative Works based upon the Original Work, and doing so without honoring these terms and 
 conditions is prohibited by copyright law and international treaty. Nothing in this License 
 is intended to affect copyright exceptions and limitations (including "fair use" or "fair 
 dealing"). This License shall terminate immediately and You may no longer exercise any of 
 the rights granted to You by this License upon your failure to honor the conditions in Section 1(c).

 10) Termination for Patent Action. This License shall terminate automatically and You 
 may no longer exercise any of the rights granted to You by this License as of the date 
 You commence an action, including a cross-claim or counterclaim, against Licensor or any 
 licensee alleging that the Original Work infringes a patent. This termination provision 
 shall not apply for an action alleging patent infringement by combinations of the Original
  Work with other software or hardware.

 11) Jurisdiction, Venue and Governing Law. Any action or suit relating to this License 
 may be brought only in the courts of a jurisdiction wherein the Licensor resides or in 
 which Licensor conducts its primary business, and under the laws of that jurisdiction 
 excluding its conflict-of-law provisions. The application of the United Nations Convention 
 on Contracts for the International Sale of Goods is expressly excluded. Any use of the Original 
 Work outside the scope of this License or after its termination shall be subject to the 
 requirements and penalties of copyright or patent law in the appropriate jurisdiction. 
 This section shall survive the termination of this License.

 12) Attorneys' Fees. In any action to enforce the terms of this License or seeking 
 damages relating thereto, the prevailing party shall be entitled to recover its costs and 
 expenses, including, without limitation, reasonable attorneys' fees and costs incurred in 
 connection with such action, including any appeal of such action. This section shall survive 
 the termination of this License.

 13) Miscellaneous. If any provision of this License is held to be unenforceable, such provision 
 shall be reformed only to the extent necessary to make it enforceable.

 14) Definition of "You" in This License. "You" throughout this License, whether in upper or 
 lower case, means an individual or a legal entity exercising rights under, and complying with 
 all of the terms of, this License. For legal entities, "You" includes any entity that controls, 
 is controlled by, or is under common control with you. For purposes of this definition, "control" 
 means (i) the power, direct or indirect, to cause the direction or management of such entity, 
 whether by contract or otherwise, or (ii) ownership of fifty percent (50%) or more of the outstanding 
 shares, or (iii) beneficial ownership of such entity.

 15) Right to Use. You may use the Original Work in all ways not otherwise restricted or conditioned
  by this License or by law, and Licensor promises not to interfere with or be responsible for such uses by You.

 16) Modification of This License. This License is Copyright © 2005 Lawrence Rosen. 
 Permission is granted to copy, distribute, or communicate this License without modification. 
 Nothing in this License permits You to modify this License as applied to the Original Work or to 
 Derivative Works. However, You may modify the text of this License and copy, distribute or communicate 
 your modified version (the "Modified License") and apply it to other original works of authorship 
 subject to the following conditions: (i) You may not indicate in any way that your Modified License 
 is the "Open Software License" or "OSL" and you may not use those names in the name of your Modified 
 License; (ii) You must replace the notice specified in the first paragraph above with the notice 
 "Licensed under <insert your license name here>" or with a notice of your own that is not confusingly 
 similar to the notice in this License; and (iii) You may not claim that your original works are open 
 source software unless your Modified License has been approved by Open Source Initiative (OSI) and 
 You comply with its license review and certification process.

 17) Non-Profit Amendment. The name of this amended version of the Open Software License ("OSL 3.0") 
 is "Non-Profit Open Software License 3.0". The original OSL 3.0 license has been amended as follows:

 (a) Licensor represents and declares that it is a not-for-profit organization that derives no revenue 
 whatsoever from the distribution of the Original Work or Derivative Works thereof, or from support 
 or services relating thereto.

 (b) The first sentence of Section 7 ["Warranty of Provenance"] of OSL 3.0 has been stricken. For 
 Original Works licensed under this Non-Profit OSL 3.0, LICENSOR OFFERS NO WARRANTIES WHATSOEVER.

 (c) In the first sentence of Section 8 ["Limitation of Liability"] of this Non-Profit OSL 3.0, 
 the list of damages for which LIABILITY IS LIMITED now includes "direct" damages.

 (d) The proviso in Section 1(c) of this License now refers to this "Non-Profit Open Software 
 License" rather than the "Open Software License". You may distribute or communicate the Original 
 Work or Derivative Works thereof under this Non-Profit OSL 3.0 license only if You make the 
 representation and declaration in paragraph (a) of this Section 17. Otherwise, You shall distribute or 
 communicate the Original Work or Derivative Works thereof only under the OSL 3.0 license and You shall 
 publish clear licensing notices so stating. Also by way of clarification, this License does not authorize 
 You to distribute or communicate works under this Non-Profit OSL 3.0 if You received them under 
 the original OSL 3.0 license.

 (e) Original Works licensed under this license shall reference "Non-Profit OSL 3.0" 
 in licensing notices to distinguish them from works licensed under the original OSL 3.0 license.
Known Bugs and Limitations:
	
	
Troubleshooting:
	Ensure GPU device being used has proper compute capability (CC).
		This program uses Unified Memory, requiring CC 3.0 or above
		(Note: only mce.cu uses the GPU device)
	Ensure the data in your input files is laid out correctly. (Refer to the docs for specifics)

Contacts:
	(Programmer) Jesse L Harder, harder-j30@webmail.uwinnipeg.ca
	
Acknowledgements:
	Supervisor: Dr. Christopher Henry, P. Eng.
	
	This research has been supported by the Natural Sciences and Engineering Research Council of Canada (NSERC) grant 418413 
	and the Tesla K40 used in this research was donated by the NVIDIA Corporation
	
	Matthew C. Schmidt, Nagiza F. Samatova, Kevin Thomas, and Byung-Hoon Park
		for their article 
	"A scalable, parallel algorithm for maximal clique enumeration" 
		which provided the base maximal clique algorithm for mce.cu.
		
	Christopher J. Henry and Sheela Ramanna 
		for their article 
	"Maximal Clique Enumeration in Finding Near Neighbourhoods"
		from which the measure used in mcMeasure.cu was obtained
	
		

