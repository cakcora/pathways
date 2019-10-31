import java.util.HashSet;
import java.util.Objects;

public class PTM {

    Gene g;
    PTMType m;
    Residue r;
    int position;
    HashSet<PTM> others = null;

    public PTM(String s) {
    //s=TUBB4B ubi K379

        String[] arr = s.split(" ");
        g = new Gene(arr[0]);
       m = new PTMType(arr[1]);
       r = new Residue(arr[2].substring(0,1));
       position = Integer.parseInt(arr[2].substring(1));
    }

    public String getName() {
        return g.getName()+" "+m.getName()+" "+r.getName()+""+position;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        PTM ptmSite = (PTM) o;
        return position == ptmSite.position &&
                Objects.equals(g, ptmSite.g) &&
                Objects.equals(m, ptmSite.m) &&
                Objects.equals(r, ptmSite.r);
    }

    @Override
    public int hashCode() {
        return Objects.hash(g, m, r, position);
    }

    public Gene getGene() {
        return g;
    }

    public void addName(PTM ptm) {
        if(this.others==null){
        this.others = new HashSet<>();
        }
        this.others.add(ptm);
    }

    public boolean hasOtherNames() {
        return this.others!=null;
    }

    public HashSet<PTM> getAllGenes() {
        return this.others;
    }
}
