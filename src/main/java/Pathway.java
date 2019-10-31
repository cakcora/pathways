import java.util.HashSet;
import java.util.Objects;
import java.util.Set;

public class Pathway {

    String name;
    Set <Gene>genes;
    public Pathway(String name) {
        this.name = name;
        genes = new HashSet();
    }

    public void addGene(Gene g) {
        genes.add(g);
    }

    public String getName() {
        return name;
    }

    public Set<Gene> getGenes() {
        return genes;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Pathway pathway = (Pathway) o;
        return Objects.equals(name, pathway.name);
    }

    @Override
    public int hashCode() {
        return Objects.hash(name);
    }

    public String getGeneNames() {
        String s=" ";
        for(Gene g:genes){
            s+=g.getName()+" ";
        }
        return s;
    }

    public boolean containsGene(Gene g) {
        return genes.contains(g);
    }

    public Set<Gene> intersect(Pathway p2) {
        Set<Gene> result = new HashSet<>(this.genes);


        result.retainAll(p2.getGenes()); //
        return result;
    }
}
