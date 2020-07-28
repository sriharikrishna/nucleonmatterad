clear all

function someval = read_plot(n,filename,yscale)
  dim=2+2*n;
  #filename='dfo_snm'
  for j=1:60
    for d=1:2
      perturbstr=int2str(j);
      if j<10
        perturbstr=sprintf('%s%s',"0",perturbstr);
      end
      filenamelocal=sprintf('%s%s%s%d%s','dfo_snm_', perturbstr,'_',d);
      % Next section is a per-nucleus summary of the different Fayans computations
      % that were run to collect all the required theoretical results.
      fid = fopen([filenamelocal,'.txt']);
      if fid >= 0
        while 1
          line = fgetl(fid);
          if line==-1
             break
          end
          b=strsplit(line,",");
    
          if length(b)==2
            j=j+1;
            probstore(j,:)=cellfun(@str2num,b(1:end));
            k=0;
          else
            k=k+1;
            datastore(k,:,j)=cellfun(@str2num,b(1:dim));
          end
        end
      end
    end
  end

  for step=[0.1 0.2]
    prob=0;
    datastore(datastore==0)=NaN;
    figure(1); clf; hold on;
    for j=1:size(probstore,1)
      if probstore(j,2)==step
        prob=prob+1;
        figure(1)
        plot(datastore(:,2,j))
        best(1,prob) = datastore(1,2,j);
        grad(1,prob) = norm(datastore(1,n+3:dim,j));
        for k=2:length(datastore(:,2,j))
            best(k,prob) = min(best(k-1,prob),datastore(k,2,j));
            grad(k,prob) = norm(datastore(k,n+3:dim,j));
        end
      end
    end

    figure(1);
    box on
    axis tight
    set(gca,'yscale',yscale)
    xlabel('iterations')
    ylabel('f value')
    title(sprintf('%s %f',filename, step))
    if step==0.1
      print('-dpng',['fig_',filename,'_','1','_1.png'])
    end
    if step==0.2
      print('-dpng',['fig_',filename,'_','2','_1.png'])
    end

    figure(2)
    plot(best)
    box on
    axis tight
    set(gca,'yscale',yscale)
    xlabel('iterations')
    ylabel('best f value')
    title(sprintf('%s %f',filename, step))
    if step==0.1
      print('-dpng',['fig_',filename,'_','1','_2.png'])
    end
    if step==0.2
      print('-dpng',['fig_',filename,'_','2','_2.png'])
    end

    figure(3)
    plot(grad)
    box on
    axis tight
    set(gca,'yscale',yscale)
    xlabel('iterations')
    ylabel('gradient norm')
    title(sprintf('%s %f',filename, step))
    if step==0.1
      print('-dpng',['fig_',filename,'_','1','_3.png'])
    end
    if step==0.2
      print('-dpng',['fig_',filename,'_','2','_3.png'])
    end
  end
  some_val = 0; % could be scalar, vector, matrix etc.
end
read_plot(2,'dfo_snm','linear')
read_plot(2,'bfg_snm','linear')
read_plot(4,'dfo_pnm','linear')
read_plot(4,'bfg_pnm','linear')
