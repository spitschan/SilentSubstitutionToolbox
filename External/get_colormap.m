function cmap=get_colormap(map,nc)

%Brewer colormaps from http://colorbrewer2.org
switch map
    case 'matlab'
        cmap=colormap('lines');
    case 'brewer1'
        if nc>9
            error('Too many color categories for brewer1 (max=9)')
        end
        cmap=[228    26    28
            55   126   184
            77   175    74
            152    78   163
            255   127     0
            255   255    51
            166    86    40
            247   129   191
            153   153   153]/255;
    case 'brewer2'
        if nc>8
            error('Too many color categories for brewer2 (max=8)')
        end
        cmap=[102	194	165
            252	141	98
            141	160	203
            231	138	195
            166	216	84
            255	217	47
            229	196	148
            179	179	179]/255;
    case 'brewer3'
        if nc>12
            error('Too many color categories for brewer3 (max=12)')
        end
        cmap=[141	211	199
            255	255	179
            190	186	218
            251	128	114
            128	177	211
            253	180	98
            179	222	105
            252	205	229
            217	217	217
            188	128	189
            204	235	197
            255	237	111]/255;
    case 'brewer_pastel'
        if nc>9
            error('Too many color categories for brewer_pastel (max=9)')
        end
        cmap=[251	180	174
            179	205	227
            204	235	197
            222	203	228
            254	217	166
            255	255	204
            229	216	189
            253	218	236
            242	242	242]/255;
    case 'brewer_dark'
        if nc>8
            error('Too many color categories for brewer1 (max=8)')
        end
        cmap=[27	158	119
            217	95	2
            117	112	179
            231	41	138
            102	166	30
            230	171	2
            166	118	29
            102	102	102]/255;
end