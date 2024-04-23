umbral=(8:62); % PTF
umbral2=(0.0000000005:0.00000000395); % delta

%Recuerda que en la matriz, el umbral va en el eje X y el umbral2 en el Y
imagesc(umbral,umbral2,matriz); %por eso aqui lo pones primero porque imagesc pide de argumentos (x,y,z)
xlabel('Media de T fagocitada por Mf (PTF) ','FontSize',12,'FontName','Arial');
ylabel('Reclutamiento (Fa1) ','FontSize',12,'FontName','Arial');
set(gca,'xaxisLocation','top')

%Aqui es para graficar los umbrales de bifurcacion y el valor nominal de
%los parametros
hold on
line([0, max(umbral)], [0.0000000005,0.0000000005], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%delta - 0.0000000005
line([0, max(umbral)], [1.68e-09,1.68e-09], 'linewidth',1.5, 'color','r'); %delta
line([0, max(umbral)], [0.00000000395,0.00000000395], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%delta + 0.00000000395
hold on
line([8,8],[0, max(umbral2)],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%PTF -
line([19.16,19.16],[0, max(umbral2)],  'linewidth',1.5, 'color','r');%PTF
line([62,62],[0, max(umbral2)],  'linewidth',1.5, 'color','k','LineStyle', "--");%PTF +

r de la derivada
disp(dfdx);