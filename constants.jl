using Interpolations
SPDLIMS = [6.7, 10.0, 11.18, 13.41, 15.65, 17.88, 20.12, 22.35, 24.59, 29.06, 31.29]
SPDS = [x*1.0 for x = 1:31]
EVS = [0.1325435315000379 0.1325435315000379 0.1325172125370564  0.13127487450241704 0.13006658393466153 0.12884587568349432 0.12753708441356393 0.12602398780994126 0.12414901473343345 0.12174774236674876 0.11873955672489289 0.11523710557911473 0.11156190105641373 0.10810381219608207 0.10513804953530827 0.10275543684880961 0.10091181822130377 0.09950690094220005 0.09843678985856666 0.09761602672072915 0.09698205530588103 0.09649286064196261 0.09612360460606814 0.09586561064086806 0.09573258071997842 0.09578778500266102 0.09620749530456038 0.09685632463286037 0.09081516940788889 0.07531010545490847 0.06775689139363002;
0.1325435315000379  0.1325435315000379  0.1325172125370564  0.13127487450241704 0.13006658393466153 0.12884587568349432 0.12753708441356393 0.12602398780994126 0.12414901473343345  0.12174774236674876 0.11873955672489289 0.11523710557911473 0.11156190105641373 0.10810381219608207 0.10513804953530827 0.10275543684880961 0.10091181822130377 0.09950690094220005 0.09843678985856666 0.09761602672072915 0.09698205530588103 0.09649286064196261 0.09612360460606814 0.09586561064086806 0.09573258071997842 0.09578778500266102 0.09620749530456038 0.09685632463286037 0.09081516940788889 0.07531010545490847 0.06775689139363002;
0.1325456464580857  0.1325456464580857  0.13251932033073394 0.13127662089083472 0.1300679206264096  0.1288466937821424  0.12753714247196457 0.12602280809214936 0.12414579380079205  0.12174141362683125 0.11872918204127808 0.11522252645280975 0.11154410382966184 0.10808452204508803 0.10511890184820835 0.10273745252817583 0.10089537894207537 0.09949196406149342 0.0984231024270075  0.09760325931355504 0.09696986453761006 0.09648090780403852 0.09611155045809097 0.09585309356637492 0.09571916656387261 0.09577284153318238 0.09619008593375646 0.09683963290307426 0.0908331710529543  0.0753248335322864  0.06775926465068466;
0.1332321722268013  0.1332321722268013  0.13320298993169505 0.13181145565012684 0.13042088283324832 0.12894418208855604 0.12721071418109808 0.12492784385011221 0.12171363003620952  0.11733816030818098 0.1121331648604402  0.10702021408376722 0.10285100067190403 0.09987674171145469 0.09788975658156569 0.09657093098037874 0.09566454553591808 0.09500426294400896 0.09449041815427167 0.09406493059670479 0.09369404690859769 0.0933580178322992  0.09304528413351508 0.09274944805591101 0.0924683233922014  0.09220578968130341 0.09198017395016037 0.09184362369853757 0.0917644072328614  0.08871424831026153 0.07335728460270675;
0.1314193498408546  0.1314193498408546  0.13139621997782278 0.13028537145139968 0.12913349915105604 0.1278004781473398  0.12599849726909504 0.12320814234260938 0.1187963043409004   0.11273919127768897 0.1063732992339558  0.10136695780754783 0.09817514302219223 0.09631748465833827 0.09522421142125151 0.09452759362514329 0.0940316648726578  0.09363822991893225 0.09329845416310448 0.09298766224934982 0.09269306024700971 0.09240781982993707 0.0921282016001703  0.09185214596226737 0.09157861960927399 0.09130745822467795 0.09103979376537594 0.09077949586904169 0.09053330083187144 0.09024381097314467 0.08858941043467777;
0.13407513041444089 0.13407513041444089 0.13402537159273148 0.13156045545361478 0.12907894069808787 0.1267328316017826  0.12458909452975955 0.12262375437587833 0.12077842710506967  0.11900190664396429 0.11726109350966354 0.11553761289450809 0.11382214090928652 0.11211022740638404 0.11039984710119563 0.10869011315455926 0.10698064456809192 0.10527126795628707 0.10356188087271953 0.10185239575287472 0.10014273538403289 0.09843286763202866 0.09672285151709809 0.09501284406515836 0.09330303091414352 0.09159347828260761 0.08988361190947651 0.08816695652843759 0.08636673788749291 0.08355601485869982 0.07316976440970993;
0.1270020777123321  0.1270020777123321  0.12705370458909734 0.12930510137216805 0.13067073733290388 0.13068094633878818 0.12936078990269934 0.12695213685921433 0.12343144977970133  0.11762506923650506 0.10749622145995702 0.09836612837594302 0.09482328159360001 0.09364513114385184 0.09308127102655817 0.09269768740820869 0.09237913830551561 0.09208720998191451 0.09180696647363398 0.09153209304399434 0.09125977573763566 0.09098871150557827 0.09071827917670977 0.09044817665492699 0.09017825645390107 0.08990845021156574 0.08963873591824173 0.08936912725082366 0.08909966039228308 0.08883018565330794 0.08855793709828139;
0.13442241987083525 0.13442241987083525 0.13448354418185682 0.13629567237719128 0.13627966127579394 0.1348585769990007  0.13243257548991194 0.12931558591015377 0.12569616497554156  0.12156021260274158 0.11603001806990604 0.10582720499285891 0.0958488284042635  0.0930461446797208  0.09234742492187202 0.09197567868009768 0.09167145807106208 0.0913879843817006  0.09111240911131112 0.09084011860231991 0.090569247637941   0.09029900342671134 0.09002904118952859 0.0897592091320106  0.08948943998719477 0.0892197041923625  0.08894998968612897 0.08868029296851358 0.08841061134351844 0.0881409054013691  0.08787080367806191;
0.13887819995151396 0.13887819995151396 0.1389036253152815  0.13917621142552147 0.13796872086618187 0.13568068750095727 0.13260642323542707 0.12894864533709574 0.1248477648551394   0.12043978111269044 0.11589659931251217 0.11110554999347974 0.10390518912668256 0.09498081967371574 0.09197705178771728 0.09131123380140183 0.09097797721667038 0.09069337857459939 0.09041893609067536 0.09014749536254232 0.08987711779593324 0.08960714221296137 0.08933732261656238 0.08906756503539968 0.08879783316509367 0.08852811281381426 0.0882583983294377  0.0879886871057796  0.08771897566394954 0.08744924780217361 0.08717940685880055;
0.1337322146407963  0.1337322146407963  0.13370479457256257 0.13196023240020283 0.129486137171202   0.12648353646289337 0.12322860496109261 0.12002718661617873 0.11713051353645793  0.1146506598032365  0.1125482524153370  0.1107028998699158  0.1089943073594413  0.1073185562828707  0.1053221250567340  0.1004048422297603  0.0919462518561366  0.0895174547359373  0.0890537870679553  0.088768348358823   0.0884970602546223  0.0882270919260370  0.0879573082970580  0.0876875678877791  0.0874178419012909  0.0871481216968185  0.0868784037115458  0.0866086853870117  0.0863389624596277  0.0860692228269227  0.0857994241922609;
0.12891195878158326 0.12891195878158326 0.12888264446726144 0.1270757289024474  0.12495293653922127 0.12277474119737045 0.12066892339059238 0.11869787871141146 0.11686436811479925  0.11513516031210463 0.113469366188368   0.11183550583059304 0.11021517537579603 0.10859945279308014 0.10698087740229961 0.10529609935029516 0.10269940440726152 0.09473297028804466 0.08897020164701387 0.08811809287352382 0.08781121008388444 0.08753894112425298 0.08726894354042625 0.08699915508340869 0.08672940833480244 0.08645967599746085 0.08618994795418819 0.08592021691836191 0.08565047221548854 0.08538068906774468 0.08511079726507009]

knots = (SPDLIMS, SPDS)
EVitp = interpolate(knots, EVS, Gridded(Linear()))
