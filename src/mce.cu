/*
 ============================================================================
 Name        : mce.cu
 Author      : Jesse Harder
 Supervisor	 : Dr. Christopher Henry, P. Eng.
 Date	     : August 23, 2016
 Version	 : 1.0
 Description : This program will find the maximal clique enumerations
 	 	 	 	 for a given Adjacency Matrix (Calculated on a GPU device)
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
 ============================================================================
 */

#include <iostream>
#include <stdlib.h>
#include <fstream>		//Read input and write output files
#include <string>		//Provides string object
#include <sstream>		//Provides methods for working with strings
#include <thrust/scan.h>//Provides parallel prefix scan algorithm

static void CheckCudaErrorAux (const char *, unsigned, const char *, cudaError_t);
#define CUDA_CHECK_RETURN(value) CheckCudaErrorAux(__FILE__,__LINE__, #value, value)

texture<int, 1, cudaReadModeElementType> data_texture;
//Number of input vectors, specified by user
unsigned int VECTORS = 0;
//Highest amount of neighbours any one vertex has, as calculated by the program
unsigned int maxNeighbours = 0;

//This class is to be inherited and will automatically cudaMallocManage (allocate in unified memory)
// when a new child class is instantiated
class ManagedStruct
{
	public:
		void *operator new(size_t len){
			//Instantiate object and return pointer to unified memory address
			void *ptr;
			cudaMallocManaged(&ptr, len);
			return ptr;
		}

		void operator delete(void *ptr){
			cudaFree(ptr);
		}
};

//This structure represents a graph node, containing the 3 arrays a graph node consists of
// When instantiated, it will automatically allocate itself to shared memory
struct graphNode : public ManagedStruct{
	//A list of ids of candidate nodes to match with
	int* cand;
	//A list of ids of nodes that would be redundant to match with
	int* cnot;
	//A list of ids in the current clique being generated
	int* compsub;
};

/******************************************************************************
 * isConnectedGPU (GPU function)
 *
 * This function will determine if one int is found in a list of other
 * ints. Each int represents the id of a vertex
 *
 * [in]:
 * 		indexA: The first vertex id, the value to be looked for
 * 		indexB: The second vertex id, the index of the list to be searched though
 * 		maxNeighbours: The length of a row in data
 *
 * [out]:
 *		None
 *
 * [return]:
 * 		True if the value is found, else false
 *
 * [notes]:
 * 		Note that this function makes use of the adjacency array which has been
 * 		bound to texture memory on the device
*******************************************************************************/
__device__ bool isConnectedGPU(int indexA, int indexB, unsigned maxNeighbours){

	// Check if the number, indexA, exists in the indexB'th row of the array, data
	unsigned index = (maxNeighbours + 1) * indexB;

	for(unsigned i = 1; i <= maxNeighbours; ++i){
		//Check if the given element in the texture array is the element we're looking for
		if(tex1Dfetch(data_texture, (index +  i)) == indexA)
			//Element was found
			return true;
	}

	//Element not found
	return false;
}

/******************************************************************************
 * createNodes (GPU kernel function)
 *
 * This function will create the next level of nodes for the MCE algorithm,
 * or output the node's compsub, if it won't produce child nodes.
 * Each thread in this algorithm represents one node
 * This function closely follows the Recursive Clique Enumeration algorithm
 * provided in
 * "A scalable, parallel algorithm for maximal clique enumeration" by
 * Matthew C. Schmidt, Nagiza F. Samatova, Kevin Thomas, and Byung-Hoon Park
 *
 * [in]:
 * 		inNode: A pointer to an array of input nodes
 * 		outNode: A pointer to an array of nodes used for output
 * 		sizeArray: A prefix summed array that tells each thread where in
 * 			outNode to write it's results
 * 		maxNeighbours: The highest number of neighbours a node can possible have
 * 		workSize: The number of nodes to be operated on
 *
 * [out]:
 *		outNode: Each node produced will be written to this array
 *
 * [return]:
 * 		Void
*******************************************************************************/
__global__ void createNodes(graphNode** inNode, graphNode** outNode, unsigned* sizeArray, unsigned maxNeighbours, unsigned workSize){

	//The id of this thread
	unsigned id = threadIdx.x + (blockDim.x * blockIdx.x);
	//fixp: The vertex in cand connected to the highest number of other vertices in cand
	int fixp;
	//cur_v: The vertex currently being worked with (initially fixp)
	int cur_v;
	//Get the location where this thread will write its output
	unsigned outputIndex;


	if(id < workSize){
		outputIndex = sizeArray[id];
		int record = -1;
		//If cand and cnot are empty, write this inNode's compsub list to the output array
		if(inNode[id]->cand[0] > -1){

			//Get the vertex in cand connected to the most other vertices in cand
			for(unsigned i = 0; i < maxNeighbours; ++i){
				int count = 0;
				for(unsigned j = 0; j < maxNeighbours; ++j){
					//Check if the current vertex i is connected to vertix j
					if(inNode[id]->cand[i] != -1 && inNode[id]->cand[j] != -1 && isConnectedGPU(inNode[id]->cand[i], inNode[id]->cand[j], maxNeighbours)){
						count++;
					}
				}
				//If the current vertex has the  most connections, then it becomes fixp
				if(count > record){
					record = count;
					fixp = inNode[id]->cand[i];
				}
			}
			//Set current vector to the vector connected to the most other vectors in cand. vector
			cur_v = fixp;

			//Track how many nodes this thread has generated so far, so we know where to write out output
			unsigned outputId = 0;

			while(cur_v != -1){

				unsigned count = 0;

				//get newnot
				//newnot = all elements in not connected to fixp
				for(unsigned i = 0; i < maxNeighbours; ++i){
					if(inNode[id]->cnot[i] > -1 && isConnectedGPU(inNode[id]->cnot[i], cur_v, maxNeighbours)){
						outNode[outputIndex + outputId]->cnot[count] = inNode[id]->cnot[i];
						count++;
					}
				}

				count = 0;
				//Get newcand
				// newCand = all nodes in cand connected to cur_v
				for(unsigned i = 0; i < maxNeighbours; ++i){
					if(inNode[id]->cand[i] > -1 && isConnectedGPU(inNode[id]->cand[i], cur_v, maxNeighbours)){
						outNode[outputIndex + outputId]->cand[count] = inNode[id]->cand[i];
						count++;
					}
				}

				//get newcompsub
				//newcompsub = oldcompsub + cur_v
				for(unsigned i = 0; i <maxNeighbours; ++i){
					if(inNode[id]->compsub[i] != -1){
						outNode[outputIndex + outputId]->compsub[i] = inNode[id]->compsub[i];
					}else{
						outNode[outputIndex + outputId]->compsub[i] = cur_v;
						break;
					}
				}

				//add cur_v to cnot
				for(unsigned i = 0; i <maxNeighbours; ++i){
					if(inNode[id]->cnot[i] == -1){
						inNode[id]->cnot[i] = cur_v;
							break;
					}
				}


				//remove cur_v from cand
				for(unsigned i = 0; i <maxNeighbours; ++i){
					if(inNode[id]->cand[i] == cur_v){
						inNode[id]->cand[i] = -1;
						break;
					}
				}

				cur_v = -1;
				//Attempt to find next cur_v (the while loop terminates if cur_v remains -1)
				//Next cur_v is the first element in cand that is not connected to fixp
				for(unsigned i = 0; i < maxNeighbours; ++i){
					if(inNode[id]->cand[i] > -1 && !isConnectedGPU(inNode[id]->cand[i], fixp, maxNeighbours)){
						outputId++;
						cur_v = inNode[id]->cand[i];
						break;
					}
				}
			}
		}
	}
	__syncthreads();

}

/******************************************************************************
 * sizeKernel (GPU kernel function)
 *
 * This function will determine how many child nodes each input node will produce.
 * Each thread handles one input node
 *
 * [in]:
 * 		nodes: A pointer to the array of input nodes
 * 		sizeArray: A pointer to the output array for the calculated sizes
 * 		maxNeighbours: The highest number of neighbours a node may have
 * 		workSize: The number of nodes in nodes to be operated on
 *
 * [out]:
 *		Each thread will write the number of children its corresponding node
 *			will produce to the thread id'th spot in sizeArray
 *
 * [return]:
 * 		Void
*******************************************************************************/
__global__ void sizeKernel(graphNode** nodes, unsigned* sizeArray, unsigned maxNeighbours, unsigned workSize){

	//Shared int array used for tracking the input node's cand vector
	extern __shared__ int shCand[];

	unsigned id = threadIdx.x + (blockDim.x * blockIdx.x);
	unsigned index = threadIdx.x  * maxNeighbours;

	//Final result, how many nodes the input node will produce
	unsigned size = 0;

	//the vertex connected to input node, connected to the most other vertices that are connected to input node
	int fixp;
	//Current node being handled
	int cur_v;
	//Store the highest amount of connections any vertex in cand has
	int record = -1;

	if(id < workSize){
		//Load cand into shared memory
		for(unsigned i = 0; i < maxNeighbours; ++i){
			shCand[index + i] = nodes[id]->cand[i];
		}

		if(shCand[index] > -1){

			//Find fixp
			//For every vertex in cand, check if that vertex is connected to any of the other vertices in cand
			for(unsigned i = 0; i < maxNeighbours; ++i){
				//Track how many matches this vertex has made
				int count = 0;

				//Check if the current vertex is in any of the other vertices in cand
				for(unsigned j = 0; j < maxNeighbours; ++j){
					if(shCand[index + i] != -1 && shCand[index + j] != -1 && isConnectedGPU(shCand[index + i], shCand[index + j], maxNeighbours)){
						count++;
					}
				}

				//If the current node has the highest count, record the count and set fixp to that node id
				if(count > record){
					record = count;
					fixp = shCand[index + i];
				}
			}


			cur_v = fixp;
			//Every iteration of this loop means one more node must be created
			while(cur_v != -1){
				size++;

				//Remove cur_v from cand
				for(unsigned i = 0; i < maxNeighbours; ++i){
					if(shCand[index + i] == cur_v)
						shCand[index + i] = -1;
				}

				cur_v = -1;
				//Find the next vertex in cand that is not connected to fixp, exit (cur_v = -1) if none are found
				for(unsigned i = 0; i < maxNeighbours; ++i){
					if(shCand[index + i] > -1 && !isConnectedGPU(shCand[index + i], fixp, maxNeighbours)){
						cur_v = shCand[index + i];
						break;
					}
				}
			}
		}
	}


	//Sync threads then write results to the output array
	__syncthreads();
	if(id < workSize){
		sizeArray[id] = size;
	}
}

/******************************************************************************
 * isConnected
 *
 * This function will determine if an integer value exists in a row of a given
 * array.
 *
 * [in]:
 * 		value: The value to be searched for
 * 		row: The row to be searched though
 * 		data: The array containing the data to be searched
 *
 * [out]:
 *		None
 *
 * [return]:
 * 		True if the value is found, else false
*******************************************************************************/
bool isConnected(unsigned value, int* data, int row){
	for(unsigned i = 1; i <= maxNeighbours; i++){
		if(data[row * (maxNeighbours + 1) + i] == value){
			return true;
		}
	}
	return false;
}

/******************************************************************************
 * writeToFile
 *
 * This function will write the compsub array values of all nodes
 * that have empty cand and cnot arrays to cliques.txt
 *
 * [in]:
 * 		data: The array of graphNodes to potentially be written out
 * 		j: The current maximum size of a compsub array
 * 		out: An output stream for writing to the text file
 * 		workSize: Number of nodes in "data"
 *
 * [out]:
 *		The compsubs of all finished nodes written to cliques.txt
 *
 * [return]:
 * 		Void
 *
 * [notes]:
 * 		If the first element of cand, cnot, or compsub are -1,
 * 			all elements are -1 in that respective array
 * 		-1's will not be written to the text file
*******************************************************************************/
void writeToFile(graphNode** data, unsigned j, std::ofstream& out, unsigned workSize){

	for(unsigned i = 0; i < workSize; ++i){
		if(data[i]->cand[0] == -1 && data[i]->cnot[0] == -1){
			for(unsigned k = 0; k < j; ++k){
				out << data[i]->compsub[k] << " ";
			}
			out << std::endl;
		}
	}
}

/******************************************************************************
 * maximalCliqueEnumeration
 *
 * This function will perform the algorithm to generate all maximal cliques.
 *
 * [in]:
 * 		structArray: Array of structs holding the first level of nodes, will
 * 			also be used to write future nodes to
 * 		initNodeCount: The initial size of struct array
 * 		intData: An array containing the initial input adjacency values
 * 		maxTileWidth: A calculated value based on shared memory limits to determine
 * 			how large out sizeKernels can be
 *
 * [out]:
 *		Calls writeToFile to write out the cliques
 *
 * [return]:
 * 		Void
*******************************************************************************/
void maximalCliqueEnumeration(graphNode** structArray, unsigned initNodeCount, int* intData, unsigned maxTileWidth, bool singleCliques){

	//Declare array for holding the counts of the amounts of nodes to be generated
	unsigned prefixArraySize = VECTORS * 3;
	unsigned* prefixSummed = new unsigned[prefixArraySize];

	//Allocate and copy node size count array to device (known as prefixSummed on host)
	unsigned* dX;
	CUDA_CHECK_RETURN(cudaMalloc((void**)&dX, prefixArraySize * sizeof(unsigned)));
	CUDA_CHECK_RETURN(cudaMemcpy(dX, prefixSummed, prefixArraySize * sizeof(unsigned), cudaMemcpyHostToDevice));

	//Allocate and copy Adjacency data to device (known as intData on host)
	int* dData;
	CUDA_CHECK_RETURN(cudaMalloc((void**)&dData, VECTORS * (maxNeighbours + 1) * sizeof(int)));
	CUDA_CHECK_RETURN(cudaMemcpy(dData, intData, VECTORS * (maxNeighbours + 1) * sizeof(int), cudaMemcpyHostToDevice));

	//As intData/dData is only ever read from the GPU (Never written to on GPU, and never read from/ written to past this point
	//	on the CPU) we bind it to texture memory on the device
	cudaBindTexture(0, data_texture, dData, VECTORS * (maxNeighbours + 1) * sizeof(int));
	printf("%s\n", cudaGetErrorString(cudaGetLastError()));

	std::cout << "Performing maximal clique enumeration ..." << std::endl;

	//Cliques are formed by a tree, generated one level at a time
	graphNode** structArray2;
	unsigned workSize = initNodeCount;

	bool swap = false;	//Determine which half of loop is run (Switches on every iteration)

	// Track number of passes through the loop
	unsigned j = 0;

	unsigned i = 0;

	//Output file stream
	std::ofstream out("cliques.txt");
	do{
		//Get launch configurations for createNodes
		unsigned dimBlock = min(workSize + 1 + (32 - ((workSize + 1) % 32)), 512);
		unsigned dimGrid = ceil((float)workSize / dimBlock);
		//Get launch configuration for sizeKernel
		unsigned sizeKernelBlock = min(maxTileWidth, dimBlock);
		unsigned sizeKernelGrid = ceil((float)workSize / sizeKernelBlock);
		unsigned sharedMemSize = maxNeighbours * sizeKernelBlock * sizeof(int);
		//std::cout << j << "th loop!\n";

		if(swap == false){
			//Launch a kernel that will determine how many nodes each of the current nodes will generate
			sizeKernel <<< sizeKernelGrid, sizeKernelBlock, sharedMemSize >>>(structArray, dX, maxNeighbours, workSize);
			//printf("%s\n", cudaGetErrorString(cudaGetLastError()));
			CUDA_CHECK_RETURN(cudaDeviceSynchronize());

			//Output the compsubs of completed nodes to a text file
			//If user does not want cliques of single nodes, skip this on the first pass through
			if(j > 0 || singleCliques == true){
				writeToFile(structArray, j + 1, out, workSize);
			}

			//Copy the node sizes back from the device
			CUDA_CHECK_RETURN(cudaMemcpy(prefixSummed, dX, (workSize + 1) * sizeof(unsigned), cudaMemcpyDeviceToHost));

			//Call thrust library to prefix sum the node counts
			thrust::exclusive_scan(prefixSummed, prefixSummed + prefixArraySize, prefixSummed);

			//Put prefix summed node sizes back on the device
			CUDA_CHECK_RETURN(cudaMemcpy(dX, prefixSummed, (workSize + 1) * sizeof(unsigned), cudaMemcpyHostToDevice));


			//Allocate the memory for the nodes to be generated
			//There should be a number of nodes equal to the sum of all the results from sizeKernel
			//If no nodes are to be generated, we are done; exit
			//std::cout << "These " << workSize << " nodes will generate " << prefixSummed[workSize] <<" more nodes!\n";
			if(prefixSummed[workSize] == 0){
				cudaFree(structArray);
				break;
			}

			structArray2 = new graphNode*[prefixSummed[workSize]];
			CUDA_CHECK_RETURN(cudaMallocManaged((void**)&(structArray2), sizeof(structArray2) * (prefixSummed[workSize])));

			for(i = 0; i < prefixSummed[workSize]; ++i){
				//Allocate node and its arrays, and initialize the arrays to -1
				structArray2[i] = new graphNode;
				structArray2[i]->cand = new int[maxNeighbours];
				CUDA_CHECK_RETURN(cudaMallocManaged((void**)&(structArray2[i]->cand), sizeof(int) * maxNeighbours));
				structArray2[i]->cnot = new int[maxNeighbours];
				CUDA_CHECK_RETURN(cudaMallocManaged((void**)&(structArray2[i]->cnot), sizeof(int) * maxNeighbours));
				structArray2[i]->compsub = new int[maxNeighbours];
				CUDA_CHECK_RETURN(cudaMallocManaged((void**)&(structArray2[i]->compsub), sizeof(int) * maxNeighbours));

				for(unsigned k = 0; k < maxNeighbours; ++k){
					structArray2[i]->cand[k] = -1;
					structArray2[i]->cnot[k] = -1;
					structArray2[i]->compsub[k] = -1;
				}
			}
			CUDA_CHECK_RETURN(cudaDeviceSynchronize());

			//Call kernel to create next level of nodes
			createNodes <<< dimGrid, dimBlock >>>(structArray, structArray2, dX, maxNeighbours, workSize);

			//Free structArray 1, so it can be reallocated for output on next pass
			cudaFree(structArray);
			//Run other half of algorithm on next pass
			swap = true;

			CUDA_CHECK_RETURN(cudaDeviceSynchronize());

		}else{
			//This segment is the same as the other half of the if statement, but structArray and structArray2 are switched

			//Launch a kernel that will determine how many nodes each of the current nodes will generate
			sizeKernel <<< sizeKernelGrid, sizeKernelBlock, sharedMemSize >>>(structArray2, dX, maxNeighbours, workSize);
			CUDA_CHECK_RETURN(cudaDeviceSynchronize());

			//Write completed nodes' compsubs to the output file
			writeToFile(structArray2, j + 1, out, workSize);

			//Copy the node sizes back from the device
			CUDA_CHECK_RETURN(cudaMemcpy(prefixSummed, dX, (workSize +1) * sizeof(unsigned), cudaMemcpyDeviceToHost));

			//Call thrust library to prefix sum the node counts
			thrust::exclusive_scan(prefixSummed, prefixSummed + prefixArraySize, prefixSummed);

			//Put prefix summed node sizes back on the device
			CUDA_CHECK_RETURN(cudaMemcpy(dX, prefixSummed, (workSize + 1) * sizeof(unsigned), cudaMemcpyHostToDevice));

			//Allocate the memory for the nodes to be generated
			//There should be a number of nodes equal to the sum of all the results from sizeKernel
			//If no nodes are to be generated, we are done; exit
			//std::cout << "These " << workSize << " nodes will generate " << prefixSummed[workSize] <<" more nodes!\n";
			if(prefixSummed[workSize] == 0){
				cudaFree(structArray2);
				break;
			}
			structArray = new graphNode*[prefixSummed[workSize]];
			CUDA_CHECK_RETURN(cudaMallocManaged((void**)&(structArray), sizeof(structArray) * (prefixSummed[workSize])));
			for(i = 0; i < prefixSummed[workSize]; ++i){

				structArray[i] = new graphNode;
				structArray[i]->cand = new int[maxNeighbours];
				CUDA_CHECK_RETURN(cudaMallocManaged((void**)&(structArray[i]->cand), sizeof(int) * maxNeighbours));
				structArray[i]->cnot = new int[maxNeighbours];
				CUDA_CHECK_RETURN(cudaMallocManaged((void**)&(structArray[i]->cnot), sizeof(int) * maxNeighbours));
				structArray[i]->compsub = new int[maxNeighbours];
				CUDA_CHECK_RETURN(cudaMallocManaged((void**)&(structArray[i]->compsub), sizeof(int) * maxNeighbours));

				for(unsigned k = 0; k < maxNeighbours; ++k){
					structArray[i]->cand[k] = -1;
					structArray[i]->cnot[k] = -1;
					structArray[i]->compsub[k] = -1;
				}
			}
			CUDA_CHECK_RETURN(cudaDeviceSynchronize());

			//Call kernel to create next level of nodes
			createNodes <<< dimGrid, dimBlock >>>(structArray2, structArray, dX, maxNeighbours, workSize);

			//Free structArray 2 so it can be used for output on next pass
			cudaFree(structArray2);

			//Run other half of algorithm on next pass
			swap = false;

			CUDA_CHECK_RETURN(cudaDeviceSynchronize());
		}
		//Set the number of nodes to be operated upon on the next pass
		workSize = prefixSummed[workSize];
		j++;

	//Loop again, unless no new nodes are to be generated (Last element of prefix sum is zero)
	}while(true);

}


/******************************************************************************
 * generateInitialNodes
 *
 * This function will generate the first level of nodes for the maximal clique
 * enumeration. It will also determine the number of nodes generated.
 *
 *
 * [in]:
 * 		structArray: A pointer to an array of graphNode structures, where the
 * 			generated nodes will be written to
 * 		intData: A pointer to the array of adjacency values
 * 		fixp: The id of the vertex/node that neighbours the most of other nodes
 *
 * [out]:
 *		structArray will have all generated nodes written to it
 *
 * [return]:
 * 		initNodeCount: The count of nodes generated
 *
 * [notes]:
 * 		This function follows the recursive clique enumeration algorithm from
 * 			"A scalable, parallel algorithm for maximal clique enumeration"
 *
*******************************************************************************/
unsigned generateInitialNodes(graphNode** structArray, int* intData, unsigned fixp){
	//cur_v is the current vertex being worked with, initially fixp
	int cur_v = fixp;
	//std::cout << "fixp: "<< fixp << std::endl;

	std::cout << "Generating initial nodes ..." << std::endl;


	//Two temporary arrays used to create the initial nodes
	int* initCand = new int[VECTORS];
	int* initNot = new int[VECTORS];
	for(unsigned i = 0; i < VECTORS; ++i){
		//Initial cand is all values 0 to VECTORS-1
		initCand[i] = i;
		//Initial not is empty (all -1)
		initNot[i] = -1;
	}


	//Pointer to array of structures used to represent graph nodes
	unsigned initNodeCount = 0;
	while(cur_v >= 0){
		//Declare new graphNode and allocate its arrays into unified memory
		structArray[initNodeCount] = new graphNode;

		structArray[initNodeCount]->cand = new int[maxNeighbours];
		CUDA_CHECK_RETURN(cudaMallocManaged((void**)&(structArray[initNodeCount]->cand), sizeof(int) * maxNeighbours));

		structArray[initNodeCount]->cnot = new int[maxNeighbours];
		CUDA_CHECK_RETURN(cudaMallocManaged((void**)&(structArray[initNodeCount]->cnot), sizeof(int) * maxNeighbours));

		structArray[initNodeCount]->compsub = new int[maxNeighbours];
		CUDA_CHECK_RETURN(cudaMallocManaged((void**)&(structArray[initNodeCount]->compsub), sizeof(int) * maxNeighbours));

		for(unsigned k = 0; k < maxNeighbours; ++k){
			//Initialize all array values to null (-1)
			structArray[initNodeCount]->cand[k] = -1;
			structArray[initNodeCount]->cnot[k] = -1;
			structArray[initNodeCount]->compsub[k] = -1;
		}

		//Fill in the 'not' values
		unsigned notCount = 0;
		for(unsigned i = 0; i < initNodeCount; ++i){
			if(initNot[i] > -1 && isConnected(initNot[i], intData, cur_v)){
				structArray[initNodeCount]->cnot[notCount] = initNot[i];
				notCount++;
			}
		}

		//Fill in the 'cand' values
		unsigned candCount = 0;
		for(unsigned i = 0; i < VECTORS; ++i){
			if(initCand[i] > -1 && isConnected(initCand[i], intData, cur_v)){
				structArray[initNodeCount]->cand[candCount] = initCand[i];
				candCount++;
			}
		}
		//Initial compsub is just cur_v
		structArray[initNodeCount]->compsub[0] = cur_v;

		//Remove cur_v from cand, add it to not
		initNot[initNodeCount] = cur_v;
		initCand[cur_v] = -1;

		//Find the next vertex to operate upon
		//This vertex must be in cand and not connected to fixp
		for(unsigned i = 0; i < VECTORS; ++i){
			if(initCand[i] > -1 && !isConnected(initCand[i], intData, fixp)){
				cur_v = i;
				break;
			}else{
				cur_v = -1;
			}
		}
		initNodeCount++;
	}

	return initNodeCount;
}

/******************************************************************************
 * genereateAdjacencyValues
 *
 * This function will populate an array with its adjacency values
 * Each row represents one vertex
 * The first element in a row is the count of how many neighbours that vertex has
 * Each following element is an int id of a neighbouring vertex
 * -1's represent nulls, or no neighbours
 * The neighbours are determines by reading in a bit matrix file, wherein
 * each 1 in a row represents a neighbour
 *
 * [in]:
 * 		intData: A pointer to the array to write output to
 * 		inFile: A string containing the name of the file to read
 *
 * [out]:
 *		intData with all adjacency values written to it
 *
 * [return]:
 * 		Void
 *
 * [notes]:
 * 		maxNeighbours is previously determined to be the most amount of
 * 		neighbours any one node can have, thus maxNeighbours + 1 is the max
 * 		size of a row in intData
 *
*******************************************************************************/
void generateAdjacencyValues(int* intData, std::string inFile){
	std::fstream inputFile(inFile.c_str(), std::ios_base::in);
	if(inputFile.fail()){
		std::cerr << "Error: Adjacency Matrix file, " << inFile.c_str() << ", could not be found." << std::endl;
		exit(1);
	}
	unsigned char c;
	unsigned i = 0;
	unsigned counter = 0;
	//Fill in the values of intData, using the input adjacency matrix file
	//Counter tracks how many neighbours this vertex has
	while(inputFile >> c){
		if(i % VECTORS == 0 && i > 0){
			intData[(maxNeighbours + 1) * ((i / VECTORS) - 1)] = counter;
			counter = 0;
		}

		// c - 48 gets numeric value of '0' and '1' chars
		if(c - 48 == 1){
			intData[(i / VECTORS) * (maxNeighbours + 1) + 1 + counter] = i % VECTORS;
			counter++;
		}
		i++;
	}

}

/******************************************************************************
 * printHelp
 *
 * Prints out all available command parameters, and a short description of each
 *
 * [in]:
 * 		None
 * [out]:
 * 		A list and description of all command parameters
 *
 * [return]:
 * 		void
 *
 *******************************************************************************/
void printHelp(){
	printf("Required Parameters:\n");
	printf("\t-v [int > 0]: Specify how many vertices the adjacency matrix has. !!Required Parameter!!\n");

	printf("Optional Parameters:\n");
	printf("\t-f [complete path and file name]: Specifies which file houses the input adjacency matrix.\n");
	printf("\t-gpu [int >= 0]: Specify which device to run GPU segments on. Requires a valid device id.\n");
	printf("\t-help: Prints out available command line parameters, then exits program\n");
}

/******************************************************************************
 * mce.cu
 *
 * This program will find all of the maximal cliques given an adjacency matrix.
 * This is an GPU based implementation of the Recursive Clique Enumerate function found in
 *  "A scalable, parallel algorithm for maximal clique enumeration" by
 *  Matthew C. Schmidt, Nagiza F. Samatova, Kevin Thomas, and Byung-Hoon Park
 *  (I recommend getting an understanding of the Recursice Clique Enumeration
 *  	algorithm presented in this article before trying to read this code)
 *
 *
 *
 * [Command line parameters]:
 * 		-v [int > 0]: The number of vertices in the input file
 * 			!!!! This parameter is mandatory at run time !!!!
 *		-gpu [int >= 0]: Specify which GPU device to run device code on
 *			(must be a valid device id)
 *		-f [full path and file name]:
 *			Allows the user to specify a file other than adjMatrix.txt
 *			to be used for input
 *		-help: Prints out all command parameters, then exits program.
 *		-sc: Specifies not to output single cliques (cliques consisting of 1 node)
 *
 * [out]:
 *		All of the maximal cliques written to
 *			cliques.txt
 *
 * [return]:
 * 		Void
 *
 * [notes]:
 * 		vertices = graph nodes
 *
*******************************************************************************/
int main(int argc, const char ** argv) {

	//Begin timing the code
	cudaEvent_t start, stop;
	float elapsedTime;
	CUDA_CHECK_RETURN(cudaEventCreate(&start));
	CUDA_CHECK_RETURN(cudaEventCreate(&stop));
	CUDA_CHECK_RETURN(cudaEventRecord(start, 0));

	//Default gpu device is 0
	unsigned device = 0;
	//Output cliques that are not connected to the graph? (Cliques of 1)
	bool singleCliques = true;
	//Default input file is adjMatrix.txt
	std::string inFile = "adjMatrix.txt";

	//Set option values for each parameter entered
	for (unsigned i = 0; i < argc; ++i) {
		if (argv[i] == std::string("-v")) {
			//Supply how many vertices are in the input file
			std::stringstream convert(argv[i + 1]);
			convert >> VECTORS;
			i++;

		}else if(argv[i] == std::string("-gpu")){
			//Specify which device to use
			std::stringstream convert(argv[i + 1]);
			convert >> device;
			i++;

		}else if(argv[i] == std::string("-f")){
			//Specify a different input file
			inFile = argv[i + 1];
			i++;

		}else if(argv[i] == std::string("-sc")){
			//Output single cliques?
			singleCliques = false;
			i++;

		}else if(argv[i] == std::string("-help")){
			//Call method to print out all available parameters, then exit program
			printHelp();
			return(0);

		}else if (i > 0) {
			std::cout << "Unknown parameter " << argv[i] << ".\n Use parameter -help for a list of available params." << std::endl;
		}
	}

	//Check for valid vector size
	if (VECTORS < 1) {
		std::cerr
				<< "The number of vectors must be > 0.\n";
		exit(1);
	}

	//Read in the data
	std::fstream inputFile(inFile.c_str(), std::ios_base::in);
	if(inputFile.fail()){
		std::cerr << "Error: Adjacency Matrix file, " << inFile.c_str() << ", could not be found." << std::endl;
		exit(1);
	}


	std::cout << "Reading from " << inFile.c_str() << " ..." << std::endl;

	//Keep track of the highest amount of neighbours any node has
	maxNeighbours = 0;

	//The vertex connected to the most other vertices
	unsigned fixp = 0;
	//Counts how many connections a given vertex has
	unsigned counter = 0;

	//Temp values and counters
	unsigned i = 0, j = 0, x = 0;
	unsigned char c;
	//Read input file to get fixp, and the number of neighbours fixp has
	while(inputFile >> c){
		x = c - 48;
		counter += x;
		//data[j] = data[j] | (x << i % 8);

		i++;
		if(i % 8 == 0 || i % VECTORS == 0){
			if(i % VECTORS == 0){
				if(counter > maxNeighbours){
					maxNeighbours = counter;
					fixp = (i - 1)/ VECTORS;
				}
				counter = 0;
			}
			j++;
		}
	}


	//Get information about the available devices
	unsigned maxTileWidth = 0;
	int nDevices;
		 cudaGetDeviceCount(&nDevices);
		  for (int i = 0; i < nDevices; i++) {
		    cudaDeviceProp prop;
		    cudaGetDeviceProperties(&prop, i);
		    //Set maxtilewidth to the max the specified device can handle
		    if(device == i){
		    	maxTileWidth = min(512.0, pow((float)2,floor(log2f(prop.sharedMemPerBlock / ((maxNeighbours) * sizeof(int))))));
		    }
		  }
	std::cout << "Using device " << device << "\n";

	//Declare an integer array for the adjacency values
	//Instead of bit string where each 1 bit represents an adjacency,
	//each vertex has a list of integer ids of other vertices
	//that are adjacent to this one. The first value of each row is
	//a count of the given vertices neighbours
	// -1s represent nulls
	//There is a row for every input vertex, and each row
	// is of length (maxNeighbours + 1) (+1 is for the count)
	int* intData = new int[VECTORS * (maxNeighbours + 1)];
	for(unsigned k = 0; k < VECTORS * (maxNeighbours + 1); ++k){
		intData[k] = -1;
	}

	//Call function to fill in intData with the adjacency values
	generateAdjacencyValues(intData, inFile);



	//Declare structure for holding the initial level of nodes, and allocate to unified memory
	graphNode** structArray = new graphNode*[VECTORS];
	CUDA_CHECK_RETURN(cudaMallocManaged((void**)&(structArray), sizeof(structArray) * VECTORS));

	//Call method to generate the first level of nodes (written to structArray)
	// and get the count of those nodes (written to initNodeCount)
	unsigned initNodeCount = generateInitialNodes(structArray, intData, fixp);

	//Call method to perform the maximal clique enumeration
	maximalCliqueEnumeration(structArray, initNodeCount, intData, maxTileWidth, singleCliques);
	std::cout << "Results written to cliques.txt\n";

	//Free memory
	free(intData);

	//Stop recording time, and print the results
	CUDA_CHECK_RETURN(cudaEventRecord(stop, 0));
	CUDA_CHECK_RETURN(cudaEventSynchronize(stop));
	CUDA_CHECK_RETURN(cudaEventElapsedTime(&elapsedTime, start, stop));
	CUDA_CHECK_RETURN(cudaEventDestroy(start));
	CUDA_CHECK_RETURN(cudaEventDestroy(stop));
	std::cout << "Elapsed time: " << elapsedTime << " ms\n";

	std::cout << "Job's done." << std::endl;
	return 0;
}

/**
 * Check the return value of the CUDA runtime API call and exit
 * the application if the call has failed.
 */
static void CheckCudaErrorAux (const char *file, unsigned line, const char *statement, cudaError_t err)
{
	if (err == cudaSuccess)
		return;
	std::cerr << statement << " returned " << cudaGetErrorString(err) << "(" << err << ") at " << file << ":" << line << std::endl;
	exit (1);
}
