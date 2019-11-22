import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class Similarity {

    public static final int CELLULAR = 1;
    HashMap<Gene,HashMap<Gene,Double>> simMap;
    HashMap<Pathway,HashMap<Pathway,Double>> pathwaySimMap;
    public double pathSimilarity(Pathway p1, Pathway p2) {

        if(pathwaySimMap.containsKey(p1)){
            if(pathwaySimMap.get(p1).containsKey(p2)){
                return pathwaySimMap.get(p1).get(p2);
            }
        }
        else {
            pathwaySimMap.put(p1,new HashMap<>());
        }
        double avg =0.0;
        int sameGene=0;
        int diffGene=0;
        for(Gene gene:p1.getGenes()){
            double thisgene=0;
            Set<Gene> genes = p2.getGenes();
            if(genes.contains(gene)){
                avg =avg+1.0;
                sameGene++;
                continue;
            }
            for(Gene gene2: genes){

                double val = getVal(gene, gene2);
                if(val>thisgene) thisgene=val;
            }
            diffGene++;
            avg = avg+thisgene;
        }
        double v = avg / p1.getGenes().size();
        pathwaySimMap.get(p1).put(p2,v);
        //System.out.println(p1.getName()+"\t"+p2.getName()+"\t"+sameGene+"\t"+diffGene+" "+v);
        return v;
    }

    private double getVal(Gene gene, Gene gene2 ) {
        double val=0.0;
        if(simMap.containsKey(gene))
        if(simMap.get(gene).containsKey(gene2)) {
            val = simMap.get(gene).get(gene2);
        }
        return val;
    }

    public void prepareSimilarity(String simfile) throws IOException {
        System.out.println(simfile);
        simMap = new HashMap<>();
        pathwaySimMap = new HashMap<>();
        HashMap<Gene, HashSet<String>> map = new HashMap<>();
        String line ="";
        BufferedReader br = new BufferedReader(new FileReader(simfile));
        while((line=br.readLine())!=null){
            line=line.replaceAll("\t\t","\t");
            String arr[] =line.split("\t");
            String region = arr[0];
            for(int i=1;i<arr.length;i++){
                Gene gene= new Gene(arr[i]);
                if(!map.containsKey(gene)){
                    map.put(gene,new HashSet<>());
                }
                map.get(gene).add(region);
            }
        }
        System.out.println(map.size()+" genes in similarity file "+simfile);

        for(Gene gene:map.keySet()){
            for(Gene gene2:map.keySet()){
                if(gene.equals(gene2)){
                    continue;
                }
                else{
                    double v = cosineSimilarity(map.get(gene),map.get(gene2));
                    if(v>0.0) {
                        if (!simMap.containsKey(gene)) {
                            simMap.put(gene, new HashMap<>());
                        }
                        if (!simMap.get(gene).containsKey(gene2)) {
                            simMap.get(gene).put(gene2, v);
                        }
                    }
                }
            }
        }
        System.out.println(simMap.size()+" genes have at least one similar gene, i.e., sim>0");
    }

    public double cosineSimilarity(HashSet<String> A, HashSet<String> B) {
        if (A == null || B == null || A.size() == 0 || B.size() == 0 || A.size() != B.size()) {
            return 0;
        }

        double sumProduct = 0.0;
        double sumASq = 0.0;
        double sumBSq = 0.0;
        HashSet<String> A2= new HashSet<String>();
        A2.addAll(A);
        sumASq = A2.size();
        A2.retainAll(B);
        sumProduct = A2.size();
        if(sumProduct==0) return 0;
        sumBSq = B.size();
        double val = 0;
        if (sumProduct == 0) {
            val= 0.0;
        }
        else val =sumProduct / (Math.sqrt(sumASq) * Math.sqrt(sumBSq));
        if(val<0||val>1.005) throw new RuntimeException("wrong sim "+val);
        return val;
    }

    public Collection<Gene> getGenes() {
    return new HashSet<Gene>(simMap.keySet());
    }
}
