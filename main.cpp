#include <iostream>
#include <fstream>
#include <string>
#include <regex> // On Linux compilers needs to be replaced with <boost:regex.h> library, as standard gcc versions have issues with this standard library
#include <vector>
#include <algorithm>
#include <sstream>
#include <iterator>


int main() {
	std::fstream FASTA("pdb_seqres.txt"); // Read in big FASTA file from local directory
	std::ofstream FASTA_Output("Results.txt"); // Create Results.txt file in local directory
	std::string DataString; // Used only in std::getline, to import lines from FASTA, later relocated to a vector
	std::vector <std::string> VectorString; // Vector to hold DataString
	std::regex ProteinID(">[a-z0-9][a-z0-9][a-z0-9][a-z0-9]_[A-Z]"); //Match >117e_A etc.
	std::regex isProteinNGlc("[A][A-Z^P][ST]|[A][A-Z][C]"); // N-Glycosylation = Asn-X-/Ser/Thr(X=anything except Proline) || Asn-X-Cys(X=Anything)
	std::regex isProteinCGlc("[W][A-Z][A-Z][W]|[W][ST][A-Z][C]"); // C-Glycosylation = Trp-X-X-Trp || Trp-Ser/Thr/X-Cys
	std::regex ProteinSignature("mol:protein"); // All proteins in this FASTA file contain this signature, all RNA/DNA structures contain "mol:na" signature

	if (FASTA.is_open() && FASTA_Output.is_open()) { // Check whether both files are open
		std::cout << "Importing the file" << std::endl; 
		while (std::getline(FASTA, DataString)) { /*std::getline is used in a while loop here, because std::getline stops when it hits blank space
			therefore useful property to seperate out protein description from protein sequence. From file_FASTA to std::string_DataString.*/
			VectorString.push_back(DataString); /*Store everything in a vector, Vector elements follow a pattern: [0] = Protein Description, [1] = Protein
				Sequence [2] = "", [3] = Another Protein Description, [4] = Another Protein Sequence, [5] = "" etc.*/
		}
		std::cout << "File imported successfully" << std::endl;
		std::cout << "Running the loop, please wait" << std::endl;
		VectorString.erase(std::remove(VectorString.begin(), VectorString.end(), ""), VectorString.end()); //Remove VectorString elements containing "" value
		for (int VectorElement = 0; VectorElement < VectorString.size(); ++VectorElement) //The main loop of the program, VectorElement is the global iterator
		{
			if (VectorElement % 2 != 0) //Odd iterator values of vector are Protein Sequences
			{
				int NGlcounter = 0; //Reset the N-Glycosylation counter at the beginning of the loop
				int CGlcounter = 0; //Reset the C-Glycosylation counter at the beginning of the loop
				FASTA_Output << "\n"; //Print out a new line in the Results.txt file to keep Protein Description from Protein Sequence seperate
				std::string TemporaryContainer1 = VectorString[VectorElement]; /*TemporaryContainer holds a copy of a single VectorString element in std::string form, 
				ensured that it indeed is protein sequence the if modulus operation above*/
				std::sregex_iterator iter1(TemporaryContainer1.begin(), TemporaryContainer1.end(), isProteinNGlc); // Establish sregex_iterator search range for isProteinNGlc 
				std::sregex_iterator end1; // Establish sregex_iterator end, necessary for the loop later one
				std::sregex_iterator iter2(TemporaryContainer1.begin(), TemporaryContainer1.end(), isProteinCGlc); // same as above, just for C-glycosylation
				std::sregex_iterator end2;
				if ((std::regex_search(VectorString[VectorElement], isProteinNGlc)) && (std::regex_search(VectorString[VectorElement], ProteinID)) == false) 
					/*if statement for "if N-Glc sequence was found, do stuff below" and an insurance check that a protein description was not passed through by accident*/
					FASTA_Output << TemporaryContainer1 << " " << std::endl;
				{
					//FASTA_Output << TemporaryContainer1 << " "; 
					for (; iter1 != end1; ++iter1) // a loop to iterate through sregex_iterator to make sure that all matches of N-Glycosylation are found
					{
						++NGlcounter; //Update N-Glycosylation counter
					}
				FASTA_Output << "Detected N-linked glycosylation (s): " << NGlcounter << std::endl;
				goto CGlc;
				}
				CGlc:if ((std::regex_search(VectorString[VectorElement], isProteinCGlc)) && (std::regex_search(VectorString[VectorElement], ProteinID)) == false)
					/*if statement for "if C-Glc sequence was found, do stuff below" and an insurance check that a protein description was not passed through by accident*/
				{
					for (; iter2 != end2; ++iter2) //a loop to iterate through sregex_iterator to make sure that all matches of C-Glycosylation are found
					{
						++CGlcounter; // Update C-Glycosylation counter 
					}
				FASTA_Output << "Detected C-linked glycosylation (s): " << CGlcounter << std::endl;
				}
				else // if no matches to C-Glc and N-Glc are not found within the sequence
				{
				FASTA_Output << "No more possible glycosylations were detected." << std::endl;
				}
			}
			if (VectorElement % 2 == 0) //Even iterator values of vector are Protein Descriptors
			{
				bool isProteinID = std::regex_search(VectorString[VectorElement], ProteinID); // Finds protein ID pattern within the protein descriptor
				bool isProtein = std::regex_search(VectorString[VectorElement], ProteinSignature); // // Finds protein signature "mol:protein" within the protein descriptor
				if (isProtein == true && isProteinID == true) // if protein description says protein and proteinID exists then
				{
					std::string TemporaryContainer2 = VectorString[VectorElement]; // create a temporary string container for a vector element of VectorString
					for (int TCi = 1; TCi < 7; ++TCi) // A loop that makes sure that only protein ID is printed in the form of: 117e_A etc
					{
						FASTA_Output << TemporaryContainer2[TCi];
					}
				TemporaryContainer2.clear(); // Clear the temporary container after the loop has finished
				}
				else // To make sure that DNA sequences that are located in the FASTA file are not processed the same way as protein sequences are
				{
					if ((isProtein == false)) /*isProtein = "mol:protein", DNA/RNA sequences have signature of "mol:na", therefore if "mol:na" is passed through here,
					this makes sure that following statement s executed*/
					{
						break; // If executed, the program flow goes back to the global 'for' loop defined in line 31.
							   //This works, because for each seperate protein from the file, protein descriptor is passed through first, if descriptor contains
							   //"mol:na" signature, then the sequence associated with the descriptor is just simply going to be skipped, and therefore not outputted to Results.txt file
					}
				}
			}
		}
		FASTA_Output.close(); //Close the results.txt file
		std::cout << "Program has finished running, please hit ENTER to exit" << std::endl; //Prompt to alert the user that the program has finished running
		std::cin.get();
	}
}