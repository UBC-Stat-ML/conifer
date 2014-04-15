# Might be better to do this in gradle,
# but the first result in google seems 
# to depend on having javacc installed,
# making the build process heavier.

# Having the code generated in the repo
# makes the build more portable; only
# those that want to regenerate jj files
# need to have javacc installed this way.

# In the future, it would be nice to have
# a gradle task that includes downloading
# javacc as part of its task. Not sure if
# it exists already.

javacc NewickParser.jj
