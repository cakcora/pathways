import java.util.logging.Logger;

public class Residue {
    String name;
    public Residue(String s) {
        if(!s.equalsIgnoreCase("S")&&!s.equalsIgnoreCase("T")&&!s.equalsIgnoreCase("Y")&&!s.equalsIgnoreCase("K")){
            Logger.getGlobal().warning("Unexpected residue type:"+
                    s);

        }
        name =s;
    }

    public String getName() {
        return name;
    }
}
