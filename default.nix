{ pkgs   ? import <nixpkgs> {},
  stdenv ? pkgs.stdenv,
  fetchurl,
  buildDocs ? false
}:

stdenv.mkDerivation rec {
  name = "pbbs-include-${version}";
  version = "v1.0";

  src = fetchurl {
    url = "https://github.com/deepsea-inria/pbbs-include/archive/${version}.tar.gz";
    sha256 = "0ikj8nffi68cd64mjzfhj2xj14ck5ps3qqf4j12676bmiadhz6hq";
  };
        
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