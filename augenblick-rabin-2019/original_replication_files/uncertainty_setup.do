/*
Goal is to discretize normal pdf
This creates two matrices:
1) gridmid is a grid of midpoints that span $span with 0 as a midpoint
2) PDF is the probability of hitting the gird with that midpoint
*/

global span=6
global num=10

global num1=$num-1
matrix grid=J($num+1,1,0)
matrix gridmid=J($num,1,0)
matrix PDF=J($num,1,0)

//grid of values and midvalues on standard normal with total of $span
forvalues i=0(1)$num {
	matrix grid[(`i'+1),1]=($span/-2)+(`i')*($span/$num)
	cap matrix gridmid[(`i'+1),1]=($span/-2)+(`i'+.5)*($span/$num)
}

//extend span to get full PDF
matrix grid[1,1]=-1000
matrix grid[($num+1),1]=+1000

//calculate PDF between midpoints
forvalues i=0(1)$num {
	cap matrix PDF[(`i'+1),1]=normal(grid[(`i'+2),1])-normal(grid[(`i'+1),1])
}


