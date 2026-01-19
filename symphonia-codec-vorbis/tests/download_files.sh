mkdir files && cd files
wget -r -N -np --reject '*index.html*' -nH --cut-dirs=3 https://people.xiph.org/~xiphmont/test-vectors/vorbis/

