"""Independent CLI oracle tests. Python 3 standard library only; not needed for demo."""
import csv
import io
import random
import subprocess
import tempfile
import unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / 'build/peptide-lens'

class PeptideTests(unittest.TestCase):
    def execute(self, normal, mutant, lengths='3', extras=()):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'n.fa').write_text(normal);(p/'m.fa').write_text(mutant)
            run=subprocess.run([str(BIN),'--normal',str(p/'n.fa'),'--mutant',str(p/'m.fa'),
                '--out',str(p/'out'),'--lengths',lengths,*extras],capture_output=True,text=True)
            files={f.name:f.read_text() for f in p.glob('out_*')}
            return run,files
    def rows(self, files, kind='candidates'):
        return list(csv.DictReader(io.StringIO(files[f'out_{kind}.csv'])))
    def test_boundary_and_short_sequences(self):
        r,f=self.execute('>s\nACDEF\n','>s\nYCDEA\n','3,9')
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual([(x['start_1based'],x['end_1based']) for x in self.rows(f)],[('1','3'),('3','5')])
    def test_no_mutation(self):
        r,f=self.execute('>s\nACDE\n','>s\nACDE\n')
        self.assertEqual(r.returncode,0);self.assertEqual(self.rows(f),[])
    def test_multi_mutation_once_per_window(self):
        r,f=self.execute('>s\nACDEF\n','>s\nAYWEF\n')
        self.assertEqual(r.returncode,0)
        rows=self.rows(f);self.assertEqual(len(rows),3)
        self.assertEqual(rows[0]['mutation_positions_1based'],'2;3')
    def test_duplicate_retains_sources_and_exact_query(self):
        r,f=self.execute('>a\nACDEF\n>b\nACDEF\n','>b\nAYDEF\n>a\nAYDEF\n',extras=('--query','AYD'))
        self.assertEqual(r.returncode,0)
        self.assertEqual(len(self.rows(f)),4)
        self.assertEqual(len(self.rows(f,'unique')),2)
        self.assertEqual(len(self.rows(f,'queries')),2)
        self.assertEqual(self.rows(f,'unique')[0]['sample_count'],'2')
    def test_sample_filter_and_csv_quotes(self):
        r,f=self.execute('>a,b\nACDEF\n>c\nACDEF\n','>c\nAYDEF\n>a,b\nAYDEF\n',extras=('--sample','a,b'))
        self.assertEqual(r.returncode,0)
        self.assertEqual({x['sample_id'] for x in self.rows(f,'queries')},{'a,b'})
    def test_lowercase_multiline_and_crlf(self):
        r,f=self.execute('>s desc\r\nac\r\ndef\r\n','>s\naydef\n','3,3')
        self.assertEqual(r.returncode,0);self.assertEqual(len(self.rows(f)),2)
    def test_reject_invalid_fasta(self):
        cases=[('', '>s\nACD\n'),('ACD','>s\nACD'),('>s\nACX','>s\nACD'),
               ('>s\nACD\n>s\nACD','>s\nACD'),('>s\n','>s\nACD'),
               ('>s\nACD','>s\nACDE'),('>s\nACD','>x\nACD'),
               ('>s\nACD','>s\nACD\n>x\nACD')]
        for n,m in cases:
            with self.subTest(n=n,m=m):
                r,_=self.execute(n,m);self.assertNotEqual(r.returncode,0)
    def test_reject_invalid_lengths_and_options(self):
        for length in ['0','-3','3,','3,,4','abc','10001']:
            with self.subTest(length=length):
                r,_=self.execute('>s\nACD','>s\nAYD',length);self.assertNotEqual(r.returncode,0)
        r,_=self.execute('>s\nACD','>s\nAYD',extras=('--sample','missing'))
        self.assertNotEqual(r.returncode,0)
    def test_randomized_independent_oracle(self):
        rng=random.Random(1729);alphabet='ACDEFGHIKLMNPQRSTVWY'
        for trial in range(80):
            n=rng.randint(1,35)
            normal=''.join(rng.choice(alphabet) for _ in range(n));mutant=list(normal)
            for i in range(n):
                if rng.random()<0.2:mutant[i]=rng.choice(alphabet)
            mutant=''.join(mutant)
            r,f=self.execute('>s\n'+normal,'>s\n'+mutant,'1,3,8')
            self.assertEqual(r.returncode,0,r.stderr)
            expected=set()
            for k in (1,3,8):
                for start in range(n-k+1):
                    a=normal[start:start+k];b=mutant[start:start+k]
                    if a!=b:
                        pos=';'.join(str(start+j+1) for j in range(k) if a[j]!=b[j])
                        expected.add((str(start+1),str(start+k),a,b,pos))
            actual={(x['start_1based'],x['end_1based'],x['normal_peptide'],x['mutant_peptide'],x['mutation_positions_1based']) for x in self.rows(f)}
            self.assertEqual(actual,expected,f'trial {trial}')
    def test_bundled_ground_truth(self):
        r,f=self.execute((ROOT/'data/normal.fasta').read_text(),(ROOT/'data/mutant.fasta').read_text(),'8,9,10,11')
        self.assertEqual(r.returncode,0)
        expected=list(csv.DictReader(io.StringIO((ROOT/'data/expected_mutations.csv').read_text())))
        self.assertEqual(self.rows(f,'mutations'),expected)
        self.assertNotIn('demo_control',{x['sample_id'] for x in self.rows(f)})

if __name__=='__main__': unittest.main(verbosity=2)
