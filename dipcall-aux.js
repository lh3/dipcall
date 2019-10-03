#!/usr/bin/env k8

var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

function vcfpair(args)
{
	var c, is_male = false, sample = 'syndip', fn_par = null, par = null;
	while ((c = getopt(args, "ms:p:")) != null) {
		if (c == 's') sample = getopt.arg;
		else if (c == 'p') fn_par = getopt.arg, is_male = true;
	}
	if (getopt.ind == args.length) {
		print("Usage: dipcall-aux.js vcfpair [options] <in.pair.vcf>");
		print("Options:");
		print("  -p FILE  chrX PAR; assuming male sample []");
		print("  -s STR   sample name [" + sample + "]");
		exit(1);
	}

	if (fn_par) {
		var file = new File(fn_par);
		var buf = new Bytes();
		par = [];
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			if (/^(chr)?X$/.test(t[0]))
				par.push([parseInt(t[1]), parseInt(t[2])]);
		}
		buf.destroy();
		file.close();
	}

	var re_ctg = is_male? /^(chr)?([0-9]+|X|Y)$/ : /^(chr)?([0-9]+|X)$/;
	var label = ['1', '2'];
	var buf = new Bytes();
	var file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var m, line = buf.toString();
		if (line.charAt(0) == '#') {
			if (/^##(source|reference)=/.test(line)) continue;
			if ((m = /^##contig=.*ID=([^\s,]+)/.exec(line)) != null) {
				if (!re_ctg.test(m[1])) continue;
			} else if (/^#CHROM/.test(line)) {
				var t = line.split("\t");
				--t.length;
				t[t.length-1] = sample;
				line = t.join("\t");
				print('##FILTER=<ID=HET1,Description="Heterozygous in the first haplotype">');
				print('##FILTER=<ID=HET2,Description="Heterozygous in the second haplotype">');
				print('##FILTER=<ID=GAP1,Description="Uncalled in the first haplotype">');
				print('##FILTER=<ID=GAP2,Description="Uncalled in the second haplotype">');
				if (is_male) {
					print('##FILTER=<ID=DIPX,Description="Diploid chrX in non-PAR">');
					print('##FILTER=<ID=DIPY,Description="Diploid chrY in non-PAR">');
				}
			}
			print(line);
			continue;
		}
		var t = line.split("\t");
		if (!re_ctg.test(t[0])) continue;
		var GT = null, AD = null, FILTER = [], HT = [null, null];
		for (var i = 0; i < 2; ++i) {
			if ((m = /^(\.|[0-9]+)\/(\.|[0-9]+):(\S+)/.exec(t[9+i])) == null) {
				warn(line);
				throw Error("malformatted VCF");
			}
			var s = m[3].split(",");
			if (AD == null) {
				AD = [];
				for (var j = 0; j < s.length; ++j)
					AD[j] = 0;
			}
			for (var j = 0; j < s.length; ++j)
				AD[j] += parseInt(s[j]);
			if (m[1] == '.') {
				FILTER.push('GAP' + label[i]);
				HT[i] = '.';
			} else if (m[1] != m[2]) {
				FILTER.push('HET' + label[i]);
				HT[i] = '.';
			} else HT[i] = m[1];
		}
		--t.length;
		// test if this is in a haploid region
		var hap = 0, st = parseInt(t[1]), en = st + t[3].length;
		if (is_male) {
			if (/^(chr)?X/.test(t[0])) {
				if (par != null) {
					var r = par, in_par = false;
					for (var i = 0; i < r.length; ++i)
						if (r[i][0] <= st && en <= r[i][1])
							in_par = true;
					hap = in_par? 0 : 2;
				}
			} else if (/^(chr)?Y/.test(t[0])) {
				hap = 1;
			}
		}
		// special treatment for haploid regions
		if (hap > 0 && FILTER.length == 1) {
			if ((hap == 2 && FILTER[0] == "GAP1") || (hap == 1 && FILTER[0] == "GAP2"))
				FILTER.length = 0;
		}
		if (hap == 2 && HT[0] != '.') FILTER.push("DIPX");
		if (hap == 1 && HT[1] != '.') FILTER.push("DIPY");
		// update VCF
		t[5] = 30; // fake QUAL
		t[6] = FILTER.length? FILTER.join(";") : ".";
		t[9] = HT.join("|") + ":" + AD.join(",");
		print(t.join("\t"));
	}
	file.close();
	buf.destroy();
}

function samflt(args)
{
	var c, min_var_len = 50000, min_mapq = 5;

	while ((c = getopt(args, "L:q:")) != null) {
		if (c == 'L') min_var_len = parseInt(getopt.arg);
		else if (c == 'q') min_mapq = parseInt(getopt.arg);
	}

	if (args.length - getopt.ind < 1) {
		warn("Usage: dipcall-aux.js samflt [-L minAlnBlockLen] [-q minMAPQ] <in.sam.gz>");
		exit(1);
	}

	var re = /(\d+)([MIDSH])/g;
	var file = new File(args[getopt.ind]);
	var buf = new Bytes();

	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		if (line.charAt(0) == '@') print(line);
		var m, t = line.split("\t", 6);
		var flag = parseInt(t[1]);
		if (flag & 0x100) continue;
		if (parseInt(t[4]) < min_mapq) continue;
		var blen = 0;
		while ((m = re.exec(t[5])) != null)
			if (m[2] == 'M' || m[2] == 'I' || m[2] == 'D')
				blen += parseInt(m[1]);
		if (blen < min_var_len) continue;
		print(line);
	}

	buf.destroy();
	file.close();
}

function main(args)
{
	if (args.length == 0) {
		print("Usage: dipcall-aux.js <command> [arguments]");
		print("Commands:");
		print("  samflt     filter SAM file");
		print("  vcfpair    convert 2-sample VCF to phased VCF");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'samflt') samflt(args);
	else if (cmd == 'vcfpair') vcfpair(args);
	else throw("unrecognized command: " + cmd);
}

main(arguments);
