#include <SourceFinder.h>
#include <ConfigParser.h>
#include <RInside.h>

#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

using namespace std;

void Usage(){
	cout<<"*** USAGE ***"<<endl;
	cout<<"[EXE] --config=[CONFIG-FILE]"<<endl;
	cout<<"*************"<<endl;
}

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "config", required_argument, 0, 'c' },
  {(char*)0, (int)0, (int*)0, (int)0}
};

std::string configFileName= "";


int main(int argc, char *argv[]){

	int c = 0;
  int option_index = 0;
	std::string inputFileName= "";

	while((c = getopt_long(argc, argv, "hc:",options_tab, &option_index)) != -1) {
    
    switch (c) {
			case 'h':
			{
      	Usage();	
				exit(0);
			}
    	case 'c':	
			{
				configFileName= std::string(optarg);	
				break;	
			}
    	default:
			{
      	Usage();	
				exit(0);
			}
    }//close switch
	}//close while
 
	//## Check config file
	if(configFileName==""){
		cerr<<"ERROR: Invalid or empty config filename, see program usage!"<<endl;
		Usage();
		exit(1);
	}

	//## Create RInside instance
	RInside* fR= new RInside;

	//## Read config file
	cout<<"INFO: Reading config from file "<<configFileName<<"..."<<endl;
	ConfigParser parser(configFileName);
	parser.ReadConfig();
	parser.Print();

	//## Source finder
	cout<<"INFO: Starting source finding ..."<<endl;
	SourceFinder finder;
	finder.Run();
	cout<<"INFO: End source finding"<<endl;
	
	return 0;

}//close main


