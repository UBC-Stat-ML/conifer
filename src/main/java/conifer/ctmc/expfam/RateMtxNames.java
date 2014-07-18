package conifer.ctmc.expfam;

public enum RateMtxNames
{
    KIMURA1980("kimura1980()"),
    ACCORDANCE("accordance()"),
    PAIR("pair()"),
    POLARITY("polarity()"),
    POLARITYSIZE("polaritySize()");

    private final String stringValue;

    private RateMtxNames(final String stringValue) {
        this.stringValue = stringValue;
    }

    public static RateMtxNames fromString(final String string) {
        for (final RateMtxNames something : values()) {
            if (something.stringValue.equals(string)) {
                return something;
            }
        }
        return null;
    }

}
