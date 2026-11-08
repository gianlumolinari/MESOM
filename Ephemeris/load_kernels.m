function load_kernels()
%LOAD_KERNELS  Load every SPICE kernel in Ephemeris/kernels/.
%
% Needed: a leap-seconds kernel (naif0012.tls), a planetary ephemeris
% (e.g. de440s.bsp), GM values (e.g. gm_de440.tpc) and radii (pck00010.tpc),
% all from https://naif.jpl.nasa.gov/pub/naif/generic_kernels/.

    here = fileparts(mfilename('fullpath'));
    files = dir(fullfile(here, 'kernels', '*.*'));
    files = files(~[files.isdir] & ~startsWith({files.name}, '.') & ~endsWith({files.name}, '.md'));
    assert(~isempty(files), 'No kernels found in %s', fullfile(here, 'kernels'));
    cspice_kclear;
    for k = 1:numel(files)
        cspice_furnsh(fullfile(files(k).folder, files(k).name));
    end
end
