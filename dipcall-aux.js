#!/usr/bin/env k8

var version = "0.1";

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
		if (/N/.test(t[4])) continue;
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

function it_index(a) {
	if (a.length == 0) return -1;
	a.sort(function(x, y) { return x[0] - y[0] });
	var last, last_i;
	for (var i = 0; i < a.length; i += 2) last = a[i][2] = a[i][1], last_i = i;
	for (var k = 1; 1<<k <= a.length; ++k) {
		var i0 = (1<<k) - 1, step = 1<<(k+1);
		for (var i = i0; i < a.length; i += step) {
			var x = 1<<(k-1);
			a[i][2] = a[i][1];
			if (a[i][2] < a[i-x][2]) a[i][2] = a[i-x][2];
			var e = i + x < a.length? a[i+x][2] : last;
			if (a[i][2] < e) a[i][2] = e;
		}
		last_i = last_i>>k&1? last_i - (1<<(k-1)) : last_i + (1<<(k-1));
		if (last_i < a.length) last = last > a[last_i][2]? last : a[last_i][2];
	}
	return k - 1;
}

function it_overlap(a, st, en) {
	var h, stack = [], b = [];
	for (h = 0; 1<<h <= a.length; ++h);
	--h;
	stack.push([(1<<h) - 1, h, 0]);
	while (stack.length) {
		var t = stack.pop();
		var x = t[0], h = t[1], w = t[2];
		if (h <= 3) {
			var i0 = x >> h << h, i1 = i0 + (1<<(h+1)) - 1;
			if (i1 >= a.length) i1 = a.length;
			for (var i = i0; i < i1 && a[i][0] < en; ++i)
				if (st < a[i][1]) b.push(a[i]);
		} else if (w == 0) { // if left child not processed
			stack.push([x, h, 1]);
			var y = x - (1<<(h-1));
			if (y >= a.length || a[y][2] > st)
				stack.push([y, h - 1, 0]);
		} else if (x < a.length && a[x][0] < en) {
			if (st < a[x][1]) b.push(a[x]);
			stack.push([x + (1<<(h-1)), h - 1, 0]);
		}
	}
	return b;
}

function dipsum(args) {
	if (args.length < 2) {
		print("Usage: dipcall-aux.js dipsum <in.dip.bed> <in.dip.vcf>");
		return;
	}
	var file, buf = new Bytes();

	var bed = {};
	file = new File(args[0]);
	var sum_len = 0;
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (bed[t[0]] == null) bed[t[0]] = [];
		var st = parseInt(t[1]), en = parseInt(t[2]);
		bed[t[0]].push([st, en]);
		sum_len += en - st; // TODO: this assumes no overlaps
	}
	file.close();
	for (var ctg in bed) it_index(bed[ctg]);

	var c = [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0]];
	file = new File(args[1]);
	while (file.readline(buf) >= 0) {
		var l = buf.toString();
		if (l[0] == '#') continue;
		var t = l.split("\t");
		if (/^(chr)?(X|Y|M)$/.test(t[0])) continue; // skip sex chromosomes
		if (t[6] != '.' && t[6] != "PASS") continue; // filtered
		if (bed[t[0]] == null) continue;
		var ref = t[3];
		var st = parseInt(t[1]) - 1;
		var en = st + ref.length;
		var a = it_overlap(bed[t[0]], st, en);
		if (a.length == 0) continue;
		// determine the variant type
		var alt = t[4].split(","), vt = null;
		for (var i = 0; i < alt.length; ++i) {
			var x = alt[i], y;
			if (x == '*') y = 3;
			if (x.length == ref.length) y = 0;
			else if (x.length > ref.length) y = 1;
			else y = 2;
			if (i == 0) vt = y;
			else if (vt != y) {
				vt = 4;
				break;
			}
		}
		// determine genotype
		var s = t[9].split(":");
		var gt = s[0].split("|");
		if (gt.length != 2) throw Error("wrong genotype");
		if (gt[0] == '.' || gt[1] == '.') continue;
		gt[0] = parseInt(gt[0]);
		gt[1] = parseInt(gt[1]);
		var k = gt[0] == gt[1]? 0 : 1;
		++c[k][vt];
	}
	file.close();
	print("Length of confident regions: " + sum_len);
	print("# Hom SNP: " + c[0][0]);
	print("# Hom INS: " + c[0][1]);
	print("# Hom DEL: " + c[0][2]);
	print("# Het SNP: " + c[1][0]);
	print("# Het INS: " + c[1][1]);
	print("# Het DEL: " + c[1][2]);
	print("# Het mixed: " + c[1][4]);
	print("SNP heterozygosity: " + (c[1][0] / sum_len).toFixed(6));
	print("Variant heterozygosity: " + ((c[1][0] + c[1][1] + c[1][2] + c[1][4]) / sum_len).toFixed(6));

	buf.destroy();
}

function main(args)
{
	if (args.length == 0) {
		print("Usage: dipcall-aux.js <command> [arguments]");
		print("Commands:");
		print("  samflt     filter SAM file");
		print("  vcfpair    convert 2-sample VCF to phased VCF");
		print("  dipsum     summarize dipcall results");
		print("  version    print version number");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'samflt') samflt(args);
	else if (cmd == 'vcfpair') vcfpair(args);
	else if (cmd == 'dipsum') dipsum(args);
	else if (cmd == 'version') print(version);
	else throw("unrecognized command: " + cmd);
}

main(arguments);
