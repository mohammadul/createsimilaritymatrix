function make()
disp('Compiling...');
p = input('Parallel?(Y/N)','s');
if(strcmpi(p,'y')==1)
    disp('Compile: P');
    mex nnfield2similarity.cpp -largeArrayDims CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -D__PARALLEL__
else
    disp('Compile: NP');
    mex nnfield2similarity.cpp -largeArrayDims
end
disp('Done.');
end
