#include <cstring>
#include <fstream>
#include <csignal>
#include <unistd.h>
using namespace std;

void signalHandler (int signum)
{
	printf("\nInterrupt signal (%d) received.\n",signum);
	exit(signum);
}

int main(int argc, char* argv[])
{
	fstream file;
	string file_name = "/dev/ttyUSB0";

	string line;
	string first_word;

	file.open(file_name,ios::out);
	if(file.is_open())
	{
		file << "execute " << argv[1] << endl;
		printf("Successfully written\n");
	}
	file.close();

	bool listening = true;
	while (listening)
	{
		if (!file.is_open())
			file.open(file_name,ios::in);
		if (getline(file,line))
		{
			first_word = line.substr(0,line.find(" "));
			line = line.substr(line.find_first_of(" \t") + 1);
			if (first_word == "executed")
			{
				printf("%s\n",line.c_str());
				while (getline(file,line))
				{
					printf("%s\n",line.c_str());
					fflush(stdout);
				}
				file.close();
				listening = false;
			}
			else
				file.close();
		}
		else
			file.close();
	}
	return 0;
}
