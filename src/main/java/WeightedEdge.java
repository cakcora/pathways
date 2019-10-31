import edu.uci.ics.jung.graph.UndirectedGraph;


public class WeightedEdge {

    int weight=0;

    public WeightedEdge(int i) {


        weight=i;
    }

    public int getWeight() {
        return weight;
    }

    public void setWeight(int weight) {
        this.weight = weight;
    }
}
