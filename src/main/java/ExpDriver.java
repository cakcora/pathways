import edu.uci.ics.jung.graph.UndirectedSparseGraph;

import java.io.FileWriter;
import java.io.IOException;

import java.util.*;
import java.util.logging.Logger;

public class ExpDriver {

    public static void main(String[] args) throws IOException {
        String dir = args[0];
        String file = "WikiPathways_2019_Human.txt";
        String clusterFile = "sites_by_cluster.txt";
        String subcellularSimFile = args[1];//GO_Cellular_Component_2018
        String molecularSimFile = args[2];//GO_Molecular_Function_2018
        String proDomainSimFile = args[3];//InterPro_Domains_2019
        String bioProcessSimFile = args[4];//GO_Biological_Process_2018
        boolean record = false;
        String d2= args[5];
        String d[] = new String[]{subcellularSimFile, molecularSimFile, proDomainSimFile, bioProcessSimFile};
        for (String fe : d) {

            int minPtmClusters = 3;
            double percOfPtms = 0.05;
            Set<Pathway> allPathways = new HashSet();
            Map<String, Cluster> clusters = new HashMap();
            Similarity sim = new Similarity();
            sim.prepareSimilarity(d2+fe);
            Set ptmsites = new HashSet<PTM>();
            Set genes = new HashSet<Gene>();

            Reader rd = new Reader();
            rd.readPathWays(allPathways, genes, dir + file);
            rd.readClusters(clusters, ptmsites, dir + clusterFile);

            HashSet<Gene> temp = new HashSet<Gene>();
            temp.addAll(genes);
            temp.retainAll(sim.getGenes());
            System.out.println(temp.size() + " genes from "+file+" are also found in the similarity file "+fe);



            int maxStep = 10;
            for (minPtmClusters = 150; minPtmClusters < 200; minPtmClusters += 20)
            //int minPtmClusters=170;
            {
                for (int step = 6; step < maxStep; step += 1)
                //int step =7;
                {
                    percOfPtms = step / (double) maxStep;
                    String fileName = dir +"networks/Cluster"+ minPtmClusters+"Perc0DOT"+step+"SimilarityIs"+fe;
                    FileWriter wr = new FileWriter(fileName+".txt");
                    UndirectedSparseGraph<Pathway, WeightedEdge> graph = new UndirectedSparseGraph();
                    for (Cluster cluster : clusters.values()) {

                        ArrayList<HashSet<Pathway>> clusterPathways = new ArrayList<>();
                        int involvedPtmCount = getInvolvedPtmCount(allPathways, cluster, clusterPathways);
                        if (involvedPtmCount / (1.0 * cluster.size()) < percOfPtms) {
                            //System.out.println("in this cluster only "+involvedPtmCount+" has connections to genes that are connected allPathways");
                            continue;
                        }
                        // System.out.println(cluster.getName()+" has "+cluster.size()+" unique ptms connects "+clusterPathways.size()+" :"+clusterPathways.toString());
                        addVerticesAndEdges(graph, clusterPathways);
                    }
                    System.out.println("With a perc threshold of "+percOfPtms+", the graph has "+graph.getVertexCount()+" vertices and "+graph.getEdgeCount()+" edges");
                    filterGraphWithThreshold(minPtmClusters, graph);
                    wr.write("pathway1\tpathway2\tsimilarity\n");
                    Pathway[] all = graph.getVertices().toArray(new Pathway[graph.getVertexCount()]);
                    double totalSim = 0.0;
                    int edge = 0;

                    for (int i = 0; i < all.length; i++) {
                        Pathway p1 = all[i];
                        for (Pathway p2 : graph.getNeighbors(p1)) {
                            edge++;
                            double simVal = sim.pathSimilarity(p1, p2);
                            if(p1.getName().compareTo(p2.getName())>0){
                                wr.write(p1.getName()+"\t"+p2.getName()+"\t"+simVal+"\r\n");
                            }
                            totalSim += simVal;
                        }
                    }
                    double avgSim = totalSim / edge;
                    if (edge == 0) avgSim = 0;
                    if (record) {
                        System.out.println(avgSim + "\t" + minPtmClusters + "\t" + percOfPtms + "\t" + graph.getVertexCount() + "\t" + graph.getEdgeCount());
                    }
                    wr.close();
                }
            }
        }
    }
    private static void filterGraphWithThreshold(int minPtmClusters, UndirectedSparseGraph<Pathway, WeightedEdge> graph) {
        int t0 = graph.getEdgeCount();
        for (WeightedEdge edge : new HashSet<WeightedEdge>(graph.getEdges())) {
            if (edge.getWeight() < minPtmClusters)
                graph.removeEdge(edge);

        }
        for (Pathway p : new HashSet<Pathway>(graph.getVertices())) {
            if (graph.getNeighborCount(p) == 0) {
                graph.removeVertex(p);
            }
            else{
                for(Gene g:p.getGenes()){
                    //System.out.println(p.getName()+"\t"+g.getName()+"\tP\tGene");
                }
                for(Pathway ne:graph.getNeighbors(p)){
                    if(p.getName().compareTo(ne.getName())>=0){
                        //System.out.println(p.getName()+"\t"+ne.getName()+"\tP\tP");
                    }
                }
            }
        }
        int t1 = graph.getEdgeCount();
        System.out.println("With a cluster threshold "+minPtmClusters+" from the graph edge count is filtered from "+t0+" to "+t1);
    }

    private static void addVerticesAndEdges(UndirectedSparseGraph<Pathway, WeightedEdge> graph, ArrayList<HashSet<Pathway>> clusterPathways) {
        for (HashSet<Pathway> h : clusterPathways) {
            for(Pathway p:h) {
                if (!graph.containsVertex(p)) {
                    graph.addVertex(p);
                }
            }
        }

        for (int i = 0; i < clusterPathways.size(); i++) {
            HashSet<Pathway> ps1 = clusterPathways.get(i);
            for (int j = i + 1; j <  clusterPathways.size(); j++) {
                HashSet<Pathway> ps2 = clusterPathways.get(j);
                for(Pathway p1:ps1){
                    for(Pathway p2:ps2){
                        WeightedEdge edge = new WeightedEdge(1);
                        if(p1.equals(p2)) continue;

                        if (graph.containsEdge(graph.findEdge(p1, p2))) {
                            edge = graph.findEdge(p1, p2);
                            edge.setWeight(1 + edge.getWeight());
                        } else graph.addEdge(edge, p1, p2);

                    }
                }
            }
        }
    }

    private static int getInvolvedPtmCount(Set<Pathway> pathways, Cluster cl, ArrayList<HashSet<Pathway>> clusterPathways) {
        int involvedPtmCount = 0;
        for (PTM ptm : cl.getPTMs()) {
            //ptm has only one gene
            HashSet<Pathway> genePathways = new HashSet<>();
            Gene gene = ptm.getGene();
            //in which pathways this gene appears?
            int count = findPathway(pathways, genePathways, gene);
            if (ptm.hasOtherNames()) {
                for (PTM p : ptm.getAllGenes()) {
                    count += findPathway(pathways, genePathways, p.getGene());
                }
            }
            if (count > 0) involvedPtmCount++;
            // else System.out.println("Gene "+gene.getName()+" is not connected to a pathway.");
            clusterPathways.add(genePathways);
        }
        if(involvedPtmCount>cl.size()) {
            Logger.getGlobal().severe(cl.size()+" involved "+involvedPtmCount);
        }
        return involvedPtmCount;
    }


    private static int findPathway(Set<Pathway> allPathways, HashSet<Pathway> pathwaysofThisGene, Gene g) {
        int count=0;
        for (Pathway p : allPathways) {
            if (p.containsGene(g)) {
                pathwaysofThisGene.add(p);
                count++;
            }
        }
        return count;
    }
}
