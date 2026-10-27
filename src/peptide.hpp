#pragma once
#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace peptide {
struct Sequence { std::string id, sequence; };
struct Mutation { std::string sample; size_t position; char reference, alternate; };
struct Candidate {
    std::string sample, mutant, normal;
    size_t start, end;
    std::vector<size_t> positions;
};
struct Analysis { std::vector<Mutation> mutations; std::vector<Candidate> candidates; };
inline std::vector<Sequence> read_fasta(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open FASTA: " + path);
    std::vector<Sequence> records;
    std::set<std::string> ids;
    std::string line;
    const std::string alphabet="ACDEFGHIKLMNPQRSTVWY";
    while (std::getline(in,line)) {
        if (!line.empty() && line.back()=='\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0]=='>') {
            std::istringstream header(line.substr(1)); std::string id; header >> id;
            if (id.empty()) throw std::runtime_error("Empty FASTA ID");
            if (!ids.insert(id).second) throw std::runtime_error("Duplicate FASTA ID: " + id);
            records.push_back({id,""});
        } else {
            if (records.empty()) throw std::runtime_error("Sequence before FASTA header");
            for (unsigned char c:line) {
                if (std::isspace(c)) continue;
                char upper=static_cast<char>(std::toupper(c));
                if (alphabet.find(upper)==std::string::npos)
                    throw std::runtime_error("Invalid amino acid in " + records.back().id + ": " + upper);
                records.back().sequence += upper;
            }
        }
    }
    if (records.empty()) throw std::runtime_error("Empty FASTA file");
    for (const auto& r:records) if (r.sequence.empty()) throw std::runtime_error("Empty sequence: " + r.id);
    return records;
}
inline Analysis analyze(const std::vector<Sequence>& normal,
                        const std::vector<Sequence>& mutant, const std::vector<size_t>& lengths) {
    if (lengths.empty()) throw std::runtime_error("No peptide lengths");
    for (size_t k:lengths) if (!k) throw std::runtime_error("Peptide length must be positive");
    std::unordered_map<std::string,std::string> lookup;
    for (const auto& r:mutant) if (!lookup.emplace(r.id,r.sequence).second)
        throw std::runtime_error("Duplicate mutant ID");
    std::set<std::string> seen;
    Analysis out;
    for (const auto& r:normal) {
        if (!seen.insert(r.id).second) throw std::runtime_error("Duplicate normal ID");
        auto it=lookup.find(r.id);
        if (it==lookup.end()) throw std::runtime_error("Missing mutant ID: " + r.id);
        const auto& m=it->second; size_t n=r.sequence.size();
        if (n!=m.size()) throw std::runtime_error("Unequal sequence lengths: " + r.id + "; substitutions only");
        std::vector<size_t> changed;
        for (size_t i=0;i<n;++i) if (r.sequence[i]!=m[i]) {
            changed.push_back(i); out.mutations.push_back({r.id,i+1,r.sequence[i],m[i]});
        }
        std::set<size_t> unique_lengths(lengths.begin(),lengths.end());
        for (size_t k:unique_lengths) {
            if (k>n) continue;
            for (size_t start=0;start<=n-k;++start) {
                auto first=std::lower_bound(changed.begin(),changed.end(),start);
                if (first==changed.end() || *first>=start+k) continue;
                Candidate c{r.id,m.substr(start,k),r.sequence.substr(start,k),start+1,start+k,{}};
                for (auto p=first;p!=changed.end() && *p<start+k;++p) c.positions.push_back(*p+1);
                out.candidates.push_back(std::move(c));
            }
        }
    }
    if (seen.size()!=lookup.size()) throw std::runtime_error("Unmatched extra mutant IDs");
    std::sort(out.candidates.begin(),out.candidates.end(),[](const Candidate& a,const Candidate& b){
        if (a.sample!=b.sample) return a.sample<b.sample;
        if (a.start!=b.start) return a.start<b.start;
        return a.end<b.end;
    });
    return out;
}
// Full occurrence records remain intact; this index only groups identical mutant peptides.
inline std::unordered_map<std::string,std::vector<size_t>> index_peptides(const Analysis& a) {
    std::unordered_map<std::string,std::vector<size_t>> index;
    for (size_t i=0;i<a.candidates.size();++i) index[a.candidates[i].mutant].push_back(i);
    return index;
}
inline std::string csv(const std::string& s) {
    std::string out="\""; for(char c:s) { if(c=='"') out+='"'; out+=c; } return out+'"';
}
inline std::string positions(const std::vector<size_t>& values) {
    std::string s; for(auto v:values) { if(!s.empty()) s+=';'; s+=std::to_string(v); } return s;
}
}
