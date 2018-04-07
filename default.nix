{ pkgs   ? import <nixpkgs> {},
  stdenv ? pkgs.stdenv,
  pbbsIncludeSrc ? ./.,
  buildDocs ? false
}:

stdenv.mkDerivation rec {
  name = "pbbs-include";

  src = pbbsIncludeSrc;
        
  installPhase = ''
    mkdir -p $out/include
    cp *.h $out/include/
  '';

  meta = {
    description = "Header files of the Problem Based Benchmark Suite.";
    license = "MIT";
    homepage = http://deepsea.inria.fr/;
  };
}