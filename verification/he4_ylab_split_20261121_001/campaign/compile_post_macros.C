void compile_post_macros(){
 if(!gSystem->CompileMacro("rootfile/lambda_he4_ylab_split_20261121_001/ExtractHe4Y.C","k"))gSystem->Exit(1);
 if(!gSystem->CompileMacro("rootfile/lambda_he4_ylab_split_20261121_001/DrawHe4Y.C","k"))gSystem->Exit(1);
 std::cout<<"HE4_Y_POST_MACROS_COMPILE_PASS"<<std::endl;
}
