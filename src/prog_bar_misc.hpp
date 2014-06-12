//As shown on:
//http://nakkaya.com/2009/11/08/command-line-progress-bar/
void printProgBar( int num_marked_points, int num_points ){
  std::string bar;

  int percent = ceil(100 * double(num_marked_points) / num_points); 

  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }

  std::cout<< "\r" "[" << bar << "] ";
  std::cout.width( 3 );
  std::cout<< percent << "%   (" << num_marked_points << "/" << num_points << ")    " << std::flush;
}
