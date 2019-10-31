import java.util.HashSet;
import java.util.Objects;

public class Cluster {
String name;
    HashSet<PTM> members;

    public Cluster(String s) {
        name = s;
        members = new HashSet<>();
    }

    public String getName() {
        return name;
    }



    public String getPTMSiteNames() {
        String s="";
        for(PTM p:members){
            s+=" "+p.getName();
        }
        return s;
    }

    public void addPTMSite(PTM ptm) {
        members.add(ptm);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Cluster cluster = (Cluster) o;
        return Objects.equals(name, cluster.name);
    }

    @Override
    public int hashCode() {
        return Objects.hash(name);
    }

    public HashSet<PTM> getPTMs() {
        return members;
    }

    public int size() {
        return members.size();
    }
}
