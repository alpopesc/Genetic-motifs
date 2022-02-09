

#ifndef TEAM_8_TESTCLASSMOTIFS_H
#define TEAM_8_TESTCLASSMOTIFS_H

class TMotifs {
public:
    TMotifs(std::string f) : file_name(f) {}
    void update(char c, int j) {
        if (buffer_sequence.size() < 7) {
            buffer_sequence.push_back(c);
        } else {
            buffer_sequence.pop_front();
            buffer_sequence.push_back(c);
        }
        counter = j;
        entire_sequence.push_back(c);
    }

    void print(){
        for (std::size_t i(0); i < buffer_sequence.size(); ++i) {
            std::cout << buffer_sequence[i];
        }
        std::cout << "\n";
    }

    std::string getSequenceName(){
        return sequence_name;
    }

    std::string getFileName(){
        return file_name;
    }

    int getCounter(){
        return counter;
    }

    std::string getEntireSequence(){
        return entire_sequence;
    }

    std::deque<char> getBuffersequences();

    void setSequenceName(const std::string& name){
        sequence_name = name;
    }



private:
    std::deque<char> buffer_sequence;
    std::vector<std::deque<char>> buffer_sequeces;
    std::string entire_sequence;
    std::string sequence_name;
    std::string file_name;
    int counter;
};

#endif 
