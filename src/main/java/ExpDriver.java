import edu.uci.ics.jung.graph.UndirectedSparseGraph;

import java.io.FileWriter;
import java.io.IOException;

import java.util.*;

public class ExpDriver {

    public static void main(String[] args) throws IOException {
        String dir = args[0];
        String paFile1 = "WikiPathways_2019_Human.txt";
        paFile1 ="bioplanet_pathway.csv";
        String clusterFile = "sites_by_cluster.txt";
        String subcellularSimFile = args[1];//GO_Cellular_Component_2018
        String molecularSimFile = args[2];//GO_Molecular_Function_2018
        String proDomainSimFile = args[3];//InterPro_Domains_2019
        String bioProcessSimFile = args[4];//GO_Biological_Process_2018

        Set<Pathway> targetFasn = loadTargets();

        boolean record = false;
        boolean rule1=false,rule2=true;
        String d2= args[5];
        String d[] = new String[]{subcellularSimFile, molecularSimFile, proDomainSimFile, bioProcessSimFile};
        for (String fe : d) {

            int minPtmClusters = 3;
            double percOfPtms = 0.05;
            double percOfPathways = 2;
            Set<Pathway> allPathways = new HashSet();
            Map<String, Cluster> clusters = new HashMap();
            Similarity sim = new Similarity();
            sim.prepareSimilarity(d2+fe);
            Set ptmsites = new HashSet<PTM>();
            Set genes = new HashSet<Gene>();


            Reader rd = new Reader();
            //rd.readPathWays(allPathways, genes, dir + paFile1);
            rd.readNIHPathWays(allPathways, genes, dir + paFile1);
            rd.readClusters(clusters, ptmsites, dir + clusterFile);

            HashSet<Gene> temp = new HashSet<Gene>();
            temp.addAll(genes);
            temp.retainAll(sim.getGenes());
            System.out.println(temp.size() + " genes from "+paFile1+" are also found in the similarity file "+fe);



            int maxStep = 10;
            for (minPtmClusters = 150; minPtmClusters < 200; minPtmClusters += 20)
            //int minPtmClusters=170;
            {
                //for (int step = 6; step < maxStep; step += 1)
                int step =7;
                {
                    percOfPtms = step / (double) maxStep;
                    String fileName = dir +"networks/Cluster"+ minPtmClusters+"Perc0DOT"+step+"SimilarityIs"+fe;
                    FileWriter wr = new FileWriter(fileName+".txt");
                    UndirectedSparseGraph<Pathway, WeightedEdge> graph = new UndirectedSparseGraph();
                    for (Cluster cluster : clusters.values()) {

                        ArrayList<HashSet<Pathway>> clusterPathways = new ArrayList<>();
                        loadClusterPathways(allPathways, cluster, clusterPathways);
                        HashSet<Pathway> inPaths = new HashSet<>();
                        for(HashSet<Pathway> s: clusterPathways){
                            for(Pathway e:s)
                                inPaths.add(e);
                        }
                        int involvedPtmCount = clusterPathways.size();
                        int involvedPathwayCount = inPaths.size();

                        //rule 1
                        if(rule1)
                        if (involvedPtmCount / (1.0 * cluster.size()) < percOfPtms) {
                            continue;
                        }
                        //rule 2
                        if(rule2)
                        if (involvedPathwayCount / (involvedPtmCount) < percOfPathways) {
                            continue;
                        }
                        // System.out.println(cluster.getName()+" has "+cluster.size()+" unique ptms connects "+clusterPathways.size()+" :"+clusterPathways.toString());
                        addVerticesAndEdges(graph, clusterPathways);
                    }
                    filterGraphWithThreshold(minPtmClusters, graph);
                    int fasnFound=0;
                    for( Pathway node:graph.getVertices()){
                        if(targetFasn.contains(node)){
                            fasnFound++;
                        }
                    }
                    System.out.println("With a perc threshold of "+percOfPtms+", the graph has "+
                            graph.getVertexCount()+" vertices and "+graph.getEdgeCount()+" edges and fasn:"+ fasnFound);

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

    private static Set<Pathway> loadTargets() {
        Set<Pathway> i = new HashSet<>();
        String r[]=new String[]{"p73 transcription factor network", "Delta Np63 pathway", "Ghrelin pathway",
                "SREBF and miR-33 in cholesterol and lipid homeostasis", "Metabolism", "Insulin signaling pathway",
                "Fatty acid biosynthesis", "AMPK signaling", "Triglyceride biosynthesis", "Fatty acyl-CoA biosynthesis",
                "Lipid and lipoprotein metabolism", "Metabolism of vitamins and cofactors",
                "Integration of energy metabolism", "Vitamin B5 (pantothenate) metabolism",
                "Fatty acid, triacylglycerol, and ketone body metabolism", "ChREBP activates metabolic gene expression"};
        for(String s:r){
            i.add(new Pathway(s));
        }
        return i;
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

    private static void loadClusterPathways(Set<Pathway> pathways, Cluster cl, ArrayList<HashSet<Pathway>> clusterPathways) {
        int involvedPtmCount = 0;
        for (PTM ptm : cl.getPTMs()) {
            //ptm has only one gene
            HashSet<Pathway> genePathways = new HashSet<>();
            Gene gene = ptm.getGene();
            //in which pathways this gene appears?
            findPathway(pathways, genePathways, gene);
            int count =0;
            if (ptm.hasOtherNames()) {
                for (PTM p : ptm.getAllGenes()) {
                    findPathway(pathways, genePathways, p.getGene());
                }
            }
            count = genePathways.size();
            if (count > 0) involvedPtmCount++;
            // else System.out.println("Gene "+gene.getName()+" is not connected to a pathway.");
            clusterPathways.add(genePathways);
        }


        return;
    }


    private static void findPathway(Set<Pathway> allPathways, HashSet<Pathway> pathwaysofThisGene, Gene g) {

        for (Pathway p : allPathways) {
            if (p.containsGene(g)) {
                pathwaysofThisGene.add(p);
               }
        }
        return;
    }
}
