#include <cstring>
#include <fstream>
#include <csignal>
#include <array>
using namespace std;

void signalHandler (int signum)
{
	printf("\nInterrupt signal (%d) received.\n",signum);

	// cleanup and close up stuff here
	// terminate program

	exit(signum);
}


int main()
{
	fstream file;
	string file_name = "test_file";

	string line;
	string previous;
	string result;
	string first_word;

   // register signal SIGINT and signal handler
   signal(SIGINT, signalHandler);

	while(true)
	{
		if(!file.is_open())
			file.open(file_name, ios::in);

		if (getline(file,line))
		{
			first_word = line.substr(0,line.find(" "));
			line = line.substr(line.find_first_of(" \t") + 1);
			file.close();
			if (line != previous && first_word == "execute")
			{

// Executes command in terminal and pipes back output
				auto pPipe = ::popen(line.c_str(),"r");
				if (pPipe == nullptr)
					throw runtime_error("Cannot open pipe");

				array<char, 256> buffer;
				result = "";

				while(!feof(pPipe))
				{
					auto bytes = fread(buffer.data(), 1, buffer.size(), pPipe);
					result.append(buffer.data(), bytes);
				}

// Puts output from pipe into file
				file.open(file_name,ios::out);
				file << "executed " << result;
				file.close();

// Makes sure not to redo the same command over and over
				previous  = line;
			}
		}
		else
			file.close();
	}


	return 0;
}
