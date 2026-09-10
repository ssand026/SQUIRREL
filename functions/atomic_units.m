function [out] = atomic_units(unit_name)
% Returns the specified SI unit in terms of Hartree atomic units.
arguments
	unit_name string {mustBeMember(unit_name, ...
		["meters","seconds","joules","eV","grams","newtons","volts","tesla"])}
end

% base units needed for Hartree
ee = 1.602176634e-19; % coulombs
hb = 1.054571818e-34; % joule-seconds
me = 9.109383702e-31; % kilograms
a0 = 5.291772109e-11; % meters
Eh = 4.359744722e-18; % joules

% other constants
cc = 2.9979245800e+08; % speed of light
aa = 7.2973525643e-03; % fine-structure constant
kB = 1.3806490000e-23; % boltzman constant

switch unit_name
	case "meters";  out = 1/a0;
	case "seconds"; out = Eh/hb;
	case "joules";  out = 1/Eh;
	case "eV";      out = ee/Eh;
	case "grams";   out = 1e-3/me;
	case "newtons"; out = a0/Eh;
	case "volts";   out = ee/Eh;
	case "tesla";   out = a0^2 * ee/hb;
end
end