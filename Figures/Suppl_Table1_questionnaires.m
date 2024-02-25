clear all;
close all;

excel = 'D:\Valeria\RSN\data\for_sharing\data_to_make_figures\KarolinskaSleepinessAnswers.xlsx';
sheet = 'Sheet1';

[sub_excel sub_excel] = xlsread(excel,sheet,'A2:A19');
KSS_eve = xlsread(excel,sheet,'B2:B19');
KSS_mor = xlsread(excel,sheet,'C2:C19');

[h p_KSS ci stats_KSS] = ttest(KSS_mor,KSS_eve)

m_eve = nanmean(KSS_eve)
m_mor = nanmean(KSS_mor)

sd_eve = nanstd(KSS_eve)
sd_mor = nanstd(KSS_mor)

%%

clear all;
close all;

excel = 'D:\Valeria\RSN\data\for_sharing\data_to_make_figures\VisualAnalogueMoodAnswers.xlsx';
sheet = 'Sheet1';

[sub_excel sub_excel] = xlsread(excel,sheet,'A3:A20');
happy_eve = xlsread(excel,sheet,'B3:B20');
sad_eve = xlsread(excel,sheet,'C3:C20');
calm_eve = xlsread(excel,sheet,'D3:D20');
tense_eve = xlsread(excel,sheet,'E3:E20');
energetic_eve = xlsread(excel,sheet,'F3:F20');
sleepy_eve = xlsread(excel,sheet,'G3:G20');

happy_mor = xlsread(excel,sheet,'H3:H20');
sad_mor = xlsread(excel,sheet,'I3:I20');
calm_mor = xlsread(excel,sheet,'J3:J20');
tense_mor = xlsread(excel,sheet,'K3:K20');
energetic_mor = xlsread(excel,sheet,'L3:L20');
sleepy_mor = xlsread(excel,sheet,'M3:M20');

happy_sad_eve = happy_eve - sad_eve;
calm_tense_eve = calm_eve - tense_eve;
energetic_sleepy_eve = energetic_eve - sleepy_eve;

happy_sad_mor = happy_mor - sad_mor;
calm_tense_mor = calm_mor - tense_mor;
energetic_sleepy_mor = energetic_mor - sleepy_mor;

[h p_happy_sad ci stats_happy_sad] = ttest(happy_sad_eve,happy_sad_mor)
[h p_calm_tense ci stats_calm_tense] = ttest(calm_tense_eve,calm_tense_mor)
[h p_energetic_sleepy ci stats_energetic_sleepy] = ttest(energetic_sleepy_eve,energetic_sleepy_mor)

m_happy_sad_eve = nanmean(happy_sad_eve)
m_happy_sad_mor = nanmean(happy_sad_mor)

m_calm_tense_eve = nanmean(calm_tense_eve)
m_calm_tense_mor = nanmean(calm_tense_mor)

m_energetic_sleepy_eve = nanmean(energetic_sleepy_eve)
m_energetic_sleepy_mor = nanmean(energetic_sleepy_mor)


sd_happy_sad_eve = nanstd(happy_sad_eve)
sd_happy_sad_mor = nanstd(happy_sad_mor)

sd_calm_tense_eve = nanstd(calm_tense_eve)
sd_calm_tense_mor = nanstd(calm_tense_mor)

sd_energetic_sleepy_eve = nanstd(energetic_sleepy_eve)
sd_energetic_sleepy_mor = nanstd(energetic_sleepy_mor)

