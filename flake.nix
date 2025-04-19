{
  description = "Flake to build + develop tfenrichr R package";

  nixConfig = {
    bash-prompt = "\[tfenrichr$(__git_ps1 \" (%s)\")\]$ ";
  };

  inputs = {
    nixpkgs.url     = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";

  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs  = import nixpkgs { inherit system; config.allowUnfree = true; };
        rpkgs = pkgs.rPackages;

        myPkg = rpkgs.buildRPackage rec {
          name    = "tfenrichr";
          pname    = "tfenrichr";
          version = "0.0.4";
          src     = ./.;

          # R‐side dependencies:
          propagatedBuildInputs = [
		  	  rpkgs.devtools
              rpkgs.BioStrings
              rpkgs.GenomicRanges
              rpkgs.SummarizedExperiment
			  rpkgs.TFMPvalue
			  rpkgs.TFBSTools
          ];

          nativeBuildInputs = [
            pkgs.R
            pkgs.pkg-config
          ];

          # C‑library dependencies
          buildInputs = [
			pkgs.bzip2
			pkgs.curl
			pkgs.gsl
			pkgs.icu75
			pkgs.libpng
			pkgs.libxml2
          ];

		  preBuild = ''
		  make build
		  #mv build $out
		  '';

          #installPhase = ''
          #  # install the tarball into $out, with dirfns present
          #  R CMD INSTALL --library=$out .
          #'';

          # enable Nix’s R-wrapper so it injects R_LD_LIBRARY_PATH
          dontUseSetLibPath = false;

          meta = with pkgs.lib; {
            description = "…";
            license     = licenses.mit;
            maintainers = [ maintainers.kewiechecki ];
          };
        };

		#myR = pkgs.rWrapper.override {
		#	packages = with rpkgs; [
		#		devtools GenomicFeatures GenomicRanges rtracklayer
		#	];
		#};

      in rec {
        # 1) allow `nix build` with no extra attr:
        defaultPackage = myPkg;

        # 2) drop you into a shell for interactive R work:
        devShells = {
          default = pkgs.mkShell {
            name = "tfenrichr-shell";
			#packages = [ myR pkgs.git ];
            buildInputs = [
              pkgs.git
              pkgs.R
              rpkgs.GenomicFeatures
              rpkgs.GenomicRanges
			  rpkgs.rtracklayer
			  pkgs.bzip2
			  pkgs.curl
              pkgs.gsl
			  pkgs.icu75
			  pkgs.libpng
			  pkgs.libxml2
            ];
            shellHook = ''
source ${pkgs.git}/share/bash-completion/completions/git-prompt.sh

export LD_LIBRARY_PATH="${pkgs.libpng.out}/lib:${pkgs.icu75.out}/lib:${pkgs.bzip2.out}/lib:${pkgs.curl.out}/lib:${pkgs.libxml2.out}/lib:${pkgs.gsl.out}/lib:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="${pkgs.libpng.out}/lib:${pkgs.icu75.out}/lib:${pkgs.bzip2.out}/lib:${pkgs.curl.out}/lib/pkgconfig:${pkgs.libxml2.out}/lib/pkgconfig:${pkgs.gsl.out}/lib/pkgconfig:$PKG_CONFIG_PATH"
            '';
          };
        };
      });
}

