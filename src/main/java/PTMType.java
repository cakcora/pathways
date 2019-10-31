import java.util.Objects;

public class PTMType {
    String name;
    public PTMType(String s) {
        name=s;
    }

    public String getName() {
        return name;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        PTMType ptmType = (PTMType) o;
        return Objects.equals(name, ptmType.name);
    }

    @Override
    public int hashCode() {
        return Objects.hash(name);
    }
}
