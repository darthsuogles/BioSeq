ReadMe

How to run the program? 
1. Make sure that your python interpreter is located in /usr/local/bin
   If yes, run
       ./align.py <match score> <mismatch score> <gap open score> <gap extend score> <sequence file 1> <sequence file 2>
   Otherwise, assume your python interpreter is somewhere in the load path
       python align.py <match score> <mismatch score> <gap open score> <gap extend score> <sequence file 1> <sequence file 2>
   Else
       This might happen if you are using windows, please go to python.org to download python, then repeat 'Otherwise' step

2. Program Output
   It might happen that there are more than one alignment (and there is really no easy way to tell which alignment you chose since we
   didn't specify how to break-tie). We print all possible alignments, with one empty line between each of them. 

3. Other algorithms
   The file contains a moderate collection of alignment algorithms, you are welcome to explore any of them. For your convenience, we 
   provided the TAGS.
