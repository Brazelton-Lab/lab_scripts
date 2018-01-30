#! /usr/bin/env python

"""
converts .data file from metawatt to 136-dimension matrix for ESOM
will process all files ending in .data
outfile is in .lrn format: http://databionic-esom.sourceforge.net/user.html#Data_files____lrn_

usage:
python metawatt-to-esom.py


Copyright:

    metawatt-to-esom  converts .data file from metawatt to 136-dimension matrix for ESOM

    Copyright (C) 2016  William Brazelton

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import glob
for filename in glob.glob('*.data'):
	print filename
	l = []
	with open(filename) as file:
		for line in file:
			if line[0] == '#': 
				cols = line.split('  ')
				if len(cols) > 136: l.append(line)
	
	# ESOM .lrn file format settings
	n = len(l)
	m = '% 137'
	s = '% 9	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1'
	var_name = '% Key	t1	t2	t3	t4	t5	t6	t7	t8	t9	t10	t11	t12	t13	t14	t15	t16	t17	t18	t19	t20	t21	t22	t23	t24	t25	t26	t27	t28	t29	t30	t31	t32	t33	t34	t35	t36	t37	t38	t39	t40	t41	t42	t43	t44	t45	t46	t47	t48	t49	t50	t51	t52	t53	t54	t55	t56	t57	t58	t59	t60	t61	t62	t63	t64	t65	t66	t67	t68	t69	t70	t71	t72	t73	t74	t75	t76	t77	t78	t79	t80	t81	t82	t83	t84	t85	t86	t87	t88	t89	t90	t91	t92	t93	t94	t95	t96	t97	t98	t99	t100	t101	t102	t103	t104	t105	t106	t107	t108	t109	t110	t111	t112	t113	t114	t115	t116	t117	t118	t119	t120	t121	t122	t123	t124	t125	t126	t127	t128	t129	t130	t131	t132	t133	t134	t135	t136'
	
	outfilename = filename + '.lrn'
	with open(outfilename,'a') as outfile:
		outfile.write('% ' + str(n) + '\n' + m + '\n' + s + '\n' + var_name + '\n')					
		for row in l: 
			row = row.replace('  ',' ')
			items = row.split(' ')
			name = items[0].split('_')
			outfile.write(name[-1]+'\t')
			for i in items[1:-1]: outfile.write(i+'\t')
			outfile.write(items[-1])
