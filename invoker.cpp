#include "invoker.hpp"
#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <sstream>
#include <sys/wait.h>
#include "cxxopts.hpp"
int main(int argc, char* argv[]) {
    MPI::Init();
    world_size = MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();
    std::string input_fifo_name = "/tmp/wtdbg-pipe";
    input_fifo_name += std::to_string(rank);
    if (rank==0) {
        
        if (mkfifo(input_fifo_name.c_str(), 0777)!=0) {
            std::cout << "Pipe creation failed on rank " << ::rank << std::endl;
            MPI::COMM_WORLD.Abort(2);
        }
        if (mkfifo("/tmp/wtdbg-alignment-pipe.gz",0777)!=0) {
            std::cout << "Pipe creation failed on rank " << ::rank << std::endl;
            MPI::COMM_WORLD.Abort(2);
        }
        
        pid_t child = fork();
        if (child==0) {
            execlp(
                "../wtdbg2", 
                "-x",
                "rs",
                "-g3.6m",
                "-t",
                "0",
                "-i",
                input_fifo_name.c_str(),
                "--load-alignments",
                "/tmp/wtdbg-alignment-pipe.gz",
                "-fo",
                "AAA",
                NULL
            );
        }
        int pipe_file = open(input_fifo_name.c_str(), O_WRONLY);
        unlink(input_fifo_name.c_str());
        if (pipe_file==-1) {
            std::cout << "Pipe open failed on rank " << ::rank << std::endl;
            MPI::COMM_WORLD.Abort(2);
        }
        int input_file = open("pacbio_filtered.fastq", O_RDONLY|O_LARGEFILE);
        if (input_file==-1) {
            std::cout << "Input file open failed on rank " << ::rank << std::endl;
            MPI::COMM_WORLD.Abort(2);
        }
        struct stat64 a;
        fstat64(input_file, &a);
        std::cout<< "rank " << rank << " ready to Bcast size";
        MPI::COMM_WORLD.Bcast(&a.st_size,1,MPI_LONG_LONG,0);
        char *buffer = new char[buffersize];
        ssize_t read_byte;

        while (read_byte=read(input_file, buffer,buffersize)) {
            write(pipe_file,buffer,read_byte);
            MPI::COMM_WORLD.Bcast(buffer, read_byte, MPI_CHAR, rank);
        }
        close(pipe_file);
        int alignment_file = open("/tmp/wtdbg-alignment-pipe.gz", O_WRONLY);
        if (alignment_file==-1) {
            std::cout << "Alignment pipe open failed on rank" << rank << std::endl;
            MPI::COMM_WORLD.Abort(2);
        }
        unlink("/tmp/wtdbg-alignment-pipe.gz");
        int tmpf = open("Ali2.gz", O_WRONLY);
        for (int i=1;i<world_size;++i) {
            std::cout<< "receiving from " << i << std::endl;
            while (true) {
                ssize_t to_read;
                MPI::Status status;
                MPI::COMM_WORLD.Recv(&to_read,1,MPI_LONG_LONG,i,0);
                if (to_read==0) break;
                MPI::COMM_WORLD.Recv(buffer,to_read,MPI_CHAR,i,1);
                write(alignment_file,buffer,to_read);
                write(tmpf, buffer,to_read);
            }
            std::cout<< "received from " << i << std::endl;
        }
        close(alignment_file);
        close (tmpf);
        int child_status;
        delete[] buffer; buffer=nullptr;
        waitpid(child, &child_status, NULL);
        MPI::COMM_WORLD.Barrier();
        
        
    } else {
        if (mkfifo(input_fifo_name.c_str(), 0777)!=0) {
            std::cout << "Input Pipe creation failed on rank " << ::rank << std::endl;
            MPI::COMM_WORLD.Abort(2);
        }
        std::string prefix = "AAA"+std::to_string(rank);

        std::string align_filename = prefix+".alignments.gz";
        if (mkfifo(align_filename.c_str(),0777)!=0) {
            std::cout << "Align Pipe creation failed on rank " << ::rank << std::endl;
            MPI::COMM_WORLD.Abort(2);
        }
        pid_t child = fork();
        if (child==0) {
            /* char *exec_argv[] = {
                "-x",
                "rs",
                "-t0",
                "-g3.6m",
                "-i",
                const_cast<char*>(input_fifo_name.c_str()),
                "--kbm-parts",
                "4",
                "--kbm-part-num",
                const_cast<char*>(std::to_string(rank-1).c_str())
            }; */
            execlp(
                "../wtdbg2-alt",
                "-x",
                "rs",
                "-g3.6m",
                "-t",
                "0",
                "-i",
                input_fifo_name.c_str(),
                "--kbm-parts",
                std::to_string(world_size-1).c_str(),
                "--kbm-part-num",
                std::to_string(rank).c_str(),
                "-fo",
                prefix.c_str(),
                NULL
            );

        }
        int pipe_file = open(input_fifo_name.c_str(), O_WRONLY);
        if (pipe_file==-1) {
            std::cout << "Pipe open failed on rank " << ::rank << std::endl;
            MPI::COMM_WORLD.Abort(2);
        }
        unlink(input_fifo_name.c_str());
        std::cout<< "rank " << rank << " ready to Bcast size";
        off64_t size;
        MPI::COMM_WORLD.Bcast(&size,1,MPI_LONG_LONG,0);
        std::cout << "rank" << rank << " received length " << size << std::endl; 
        char *buffer = new char[buffersize];
        ssize_t read_byte=0;
        int tmp2 = open(("IN"+std::to_string(rank)).c_str(),O_CREAT|O_WRONLY);
        while (read_byte<size) {
            ssize_t to_read = size-read_byte>buffersize?buffersize:size-read_byte;
            MPI::COMM_WORLD.Bcast(buffer,to_read,MPI_CHAR,0);
            write(pipe_file, buffer, to_read);
            write(tmp2, buffer, to_read);
            read_byte+=to_read;
        }
        close(pipe_file);
        int align_file = open(align_filename.c_str(), O_RDONLY);
        if (pipe_file==-1) {
            std::cout << "Pipe open failed on rank " << ::rank << std::endl;
            MPI::COMM_WORLD.Abort(2);
        }
        unlink(align_filename.c_str());
        while (read_byte = read(align_file,buffer,buffersize)) {
            MPI::COMM_WORLD.Send(&read_byte,1,MPI_LONG_LONG,0,0);
            MPI::COMM_WORLD.Send(buffer,read_byte,MPI_CHAR,0,1);

        }
        MPI::COMM_WORLD.Send(&read_byte,1,MPI_LONG_LONG,0,0);
        delete[] buffer;buffer=nullptr;
        int child_status;
        waitpid(child, &child_status, NULL);
        MPI::COMM_WORLD.Barrier();
        
    }
    MPI::Finalize();
    return 0;
}