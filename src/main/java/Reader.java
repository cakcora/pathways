import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class Reader {

    public void readPathWays(Set<Pathway> pathways, Set genes, String s) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(s));
        String line="";
        while((line=br.readLine())!=null){

            String arr[] = line.replaceAll("\t\t","\t").split("\t");
            int l = arr.length;
            Pathway p= new Pathway(arr[0]);

            for(int i=1;i<l;i++){

                Gene g = new Gene(arr[i]);
                if(!genes.contains(g)) genes.add(g);
                p.addGene(g);
            }
            if(!pathways.contains(p))pathways.add(p);

        }
    }

    public void readNIHPathWays(Set<Pathway> pathways, Set genes, String paFile) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(paFile));
        String line="";
        br.readLine();
        HashMap<String, Pathway> ps = new HashMap();
        while((line=br.readLine())!=null){

            String arr[] = line.replaceAll("\"","").split(",");
            String pathwayName = arr[1];
            String gene = arr[3];
            Pathway pathway = null;
            if(ps.containsKey(pathwayName))
                pathway = ps.get(pathwayName);
            else {
                pathway= new Pathway(pathwayName);

            }
            Gene g = new Gene(gene);
            if(!genes.contains(g)) genes.add(g);
            pathway.addGene(g);
            ps.put(pathwayName,pathway);
        }
        pathways.addAll(ps.values());
    }

    public void readClusters(Map<String,Cluster> clusters, Set ptmsites, String file) throws IOException {
        String line="";
        BufferedReader br  = new BufferedReader(new FileReader(file));
        br.readLine();//header
        while((line=br.readLine())!=null){
            String s = line.replaceAll("\t\t", "\t");
            s = s.replaceAll("  "," ");
            String arr[] = s.split("\t");
            String cname = arr[1];
            Cluster cluster;
            if(clusters.containsKey(cname)){
                cluster = clusters.get(cname);
            }
            else{
                cluster= new Cluster(cname);
                clusters.put(cname,cluster);
            }

            String[] arr2= arr[0].split(";");
            PTM ptm = new PTM(arr2[0].trim());
            if(arr2.length==1){
                // nothing. Ptm has one name only
            }
            else{

                for(int i=1;i<arr2.length;i++){
                    PTM ptmOther = new PTM(arr2[i].trim());
                    ptm.addName(ptm);
                }
            }
            if(!ptmsites.contains(ptm)) ptmsites.add(ptm);
            cluster.addPTMSite(ptm);
        }
    }
}
