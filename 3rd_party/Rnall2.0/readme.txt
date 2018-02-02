README for Rnall

Rnall (RNA Local secondary structure prediction by Local symmetric mapping) version 2.0 is a new version RNA Local secondary structure prediction tool. Rnall 2.0 predicts local RNA secondary structures for multiple sequences in the Fasta format or scans local secondary structure of a single
genomic sequence in the Fasta format. 

Rnall scans the RNA sequence with a sliding window and extracts local secondary structure candidates based on dynamic programming. Different from Rnall 1.1, Rnall 2.0 extracts all LSS with sizes no longer than a window size. Furthermore, Rnall 2.0 incorporates thermodynamic parameters from the nearest neighbor model as the score function. Rnall has a time complexity of O(W3L) in the worse case, where W is the window size (typically a small constant less than 100 nt) and L is the sequence length. In practice, we observed an average computing time in O(W2L).

Rnall has various potential applications, such as local secondary structure prediction in RNA molecules and RNA motif prediction (such as rho-independent/intrinsic terminator, riboswitch, siRNA, and viral RNA motifs). 

==============================================================================

   Rnall Version 2.0 was jointly developed by Wan, X.-F., Lin, G., and Xu, D. 
   The license belongs to Digital Biology Laboratory,
   University of Missouri-Columbia. However, the software Rnall is freely 
   available for academic or non-commercial purpose at 
   http://digbio.missouri.edu/~wanx/Rnall. 

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
   FITNESS FOR A PARTICULAR PURPOSE.

   Any bug should refer to Dong Xu by xudong@missouri.edu.
   
==============================================================================
The current available version is Windows or Linux version.

---Windows
Use any unzip software to extract the files

---Linux
unzip the file:
	gunzip Rnall.tar.Z
	tar -xvf Rnall.tar

==============================================================================
Running command:

---Windows
	Rnall2.0.exe -i filename -o outputfilename -g/-m [-w 50] [-l 10] [-e -5] -c
---Linux
./Rnall2.0.exe -i filename -o outputfilename -g/-m [-w 50] [-l 10] [-e -5] -c

	-i input file. the input file should be in the Fasta format
	-o output file for predicted RNA local secondary structure
	-g for scanning a single genome in the Fasta format
	-m for multiple sequences in the Fasta format
	-w for sliding window size (default 30)
	-l for the lower window boundary for RNA LSS output (default =10, e.g. Rnall will output all of RNA LSS with size between 10 and 30 as default)
	-e for energy threshold in the unit of Kcal/mol (default = -5)

==============================================================================
Example1 for predicting RNA LSS energy landscape
test.fasta
Windows command: Rnall2.0.exe -i test.fasta -o test.out -w 50 –l 30 –e 100 -g
Linux command: ./Rnall2.0.exe -i test.fasta -o test.out -w 50 –l 30 –e 100 -g
output: test.out

Example2 for outputting RNA LSS in CT format
test.fasta
Windows command: Rnall2.0.exe -i test.fasta -o test.ct -w 50 –e -10 –g -c
Linux command: ./Rnall2.0.exe -i test.fasta -o test.ct -w 50  –e -10 –g -c
output: test.ct

