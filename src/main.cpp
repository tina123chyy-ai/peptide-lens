#include "peptide.hpp"
#include <iostream>
#include <limits>
using namespace peptide;
void help() {
    std::cout << "Peptide Lens 1.0 | Synthetic-data educational tool\n"
    "Usage: peptide-lens --normal FILE --mutant FILE --out PREFIX [--lengths 8,9,10,11]\n"
    "                    [--sample ID] [--query EXACT_PEPTIDE]\n"
    "Output directory must already exist. Coordinates are 1-based, inclusive.\n"
    "Filters affect queries.csv only; main analysis outputs always contain all samples.\n";
}
int main(int argc,char** argv) {
    try {
        if(argc==1) { help(); return 0; }
        std::map<std::string,std::string> opts;
        const std::set<std::string> allowed={"--normal","--mutant","--out","--lengths","--sample","--query"};
        for(int i=1;i<argc;++i) {
            std::string key=argv[i]; if(key=="--help") {help();return 0;}
            if(!allowed.count(key)) throw std::runtime_error("Unknown option: "+key);
            if(i+1>=argc) throw std::runtime_error("Missing value: "+key);
            if(!opts.emplace(key,argv[++i]).second) throw std::runtime_error("Repeated option: "+key);
        }
        for(const auto& key:{"--normal","--mutant","--out"})
            if(!opts.count(key) || opts[key].empty()) throw std::runtime_error(std::string("Required: ")+key);
        std::string lengths=opts.count("--lengths")?opts["--lengths"]:"8,9,10,11";
        if(lengths.empty() || lengths.back()==',') throw std::runtime_error("Invalid lengths");
        std::istringstream input(lengths); std::string part; std::vector<size_t> ks;
        while(std::getline(input,part,',')) {
            if(part.empty() || part.find_first_not_of("0123456789")!=std::string::npos)
                throw std::runtime_error("Lengths must be comma-separated positive integers");
            auto k=std::stoull(part);
            if(k==0 || k>10000) throw std::runtime_error("Length must be between 1 and 10000");
            ks.push_back(static_cast<size_t>(k));
        }
        auto normal=read_fasta(opts["--normal"]), mutant=read_fasta(opts["--mutant"]);
        auto result=analyze(normal,mutant,ks); auto index=index_peptides(result);
        if(opts.count("--sample") && std::none_of(normal.begin(),normal.end(),[&](const Sequence& s){return s.id==opts["--sample"];}))
            throw std::runtime_error("Unknown sample ID");
        std::string prefix=opts["--out"];
        std::vector<std::ofstream> files;
        for(const auto& suffix:{"_mutations.csv","_candidates.csv","_unique.csv","_queries.csv","_summary.txt"}) {
            files.emplace_back(prefix+suffix);
            if(!files.back()) throw std::runtime_error("Cannot write output: "+prefix+suffix);
        }
        auto& mutations=files[0]; auto& candidates=files[1]; auto& unique=files[2]; auto& queries=files[3]; auto& summary=files[4];
        mutations<<"sample_id,position_1based,normal_aa,mutant_aa\n";
        for(const auto& m:result.mutations) mutations<<csv(m.sample)<<','<<m.position<<','<<m.reference<<','<<m.alternate<<'\n';
        const std::string header="sample_id,start_1based,end_1based,length,normal_peptide,mutant_peptide,mutation_positions_1based\n";
        candidates<<header; queries<<header;
        auto row=[](std::ostream& o,const Candidate& c) {
            o<<csv(c.sample)<<','<<c.start<<','<<c.end<<','<<c.mutant.size()<<','<<c.normal<<','<<c.mutant<<','<<csv(positions(c.positions))<<'\n';
        };
        for(const auto& c:result.candidates) row(candidates,c);
        size_t hits=0;
        auto emit=[&](const Candidate& c) { if(!opts.count("--sample") || c.sample==opts["--sample"]) {row(queries,c);++hits;} };
        if(opts.count("--query")) {
            auto found=index.find(opts["--query"]);
            if(found!=index.end()) for(size_t i:found->second) emit(result.candidates[i]);
        } else if(opts.count("--sample")) {for(const auto& c:result.candidates) emit(c);}
        unique<<"mutant_peptide,length,occurrence_count,sample_count\n";
        std::vector<std::string> keys; for(const auto& kv:index) keys.push_back(kv.first); std::sort(keys.begin(),keys.end());
        for(const auto& key:keys) {
            std::set<std::string> ids; for(auto i:index.at(key)) ids.insert(result.candidates[i].sample);
            unique<<key<<','<<key.size()<<','<<index.at(key).size()<<','<<ids.size()<<'\n';
        }
        std::ostringstream report;
        report<<"Peptide Lens 1.0\nSamples: "<<normal.size()<<"\nSubstitutions: "<<result.mutations.size()
        <<"\nCandidate occurrences: "<<result.candidates.size()<<"\nUnique mutant peptides: "<<index.size()
        <<"\nQuery hits: "<<hits<<"\nRequested lengths: "<<lengths
        <<"\nCoordinates: 1-based inclusive\nNo HLA binding or immunogenicity prediction is performed.\n";
        summary<<report.str();
        for(auto& f:files) {f.flush();if(!f) throw std::runtime_error("Output write failed");}
        std::cout<<report.str()<<"Saved outputs with prefix: "<<prefix<<'\n';
        return 0;
    } catch(const std::exception& e) {std::cerr<<"Error: "<<e.what()<<'\n';return 1;}
}
