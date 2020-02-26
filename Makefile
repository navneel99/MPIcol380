compile:
	mpicxx -std=c++11 -o b_m.out block_matmul.cpp
	mpicxx -std=c++11 -o nb_m.out nonblock_matmul.cpp
	mpicxx -std=c++11 -o col_m.out col_matmul.cpp

clean:
	rm *.out
	