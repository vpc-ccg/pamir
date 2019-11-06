#include "process_range.h"
#include "process_orphans.h"
#include <iostream>
#include <string>
int main(int argc, char **argv){
    
    if (argc == 1){
        std::cerr << "./process [range, orphan] [args]" << std::endl;
        return -1;
    }
    std::string tool(argv[1]);
    if(tool == "range"){
        return cram_process::main(argc-1,argv+1);
    }
    if(tool == "orphan"){
        return orphan_process::main(argc-1,argv+1);
    }
    return -1;
}
