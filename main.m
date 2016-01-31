import com.mongodb.BasicDBObject;
import java.util.regex.Pattern;
import com.mongodb.util.JSON;
import com.mongodb.MongoClient;
mongoClient = MongoClient('10.91.53.225',27017);
db = mongoClient.getDB('LINCS_L1000');
readColl = db.getCollection('meta2014');
writeColl = db.getCollection('cpcd-v1.2');
sigDB = mongoClient.getDB('LINCS_L1000_limma');
sigColl = sigDB.getCollection('siginfo');

filter = BasicDBObject();
% filter only CPC batch
filter.append('det_plate',Pattern.compile('CPD'));
batches = readColl.distinct('batch',filter);
batches = j2m(batches);

bigMatPath = 'D:\Qiaonan backup\LINCS data\newData\q2norm_n1328098x22268.gctx';
load('D:\Qiaonan backup\LINCS data\newData\id2gene');
load('D:\Qiaonan backup\LINCS data\newData\newMatRid');
geneSymbols = cell(22268,1);
lmIdx = false(22268,1);
for i = 1:numel(rid)
    geneSymbols{i} = dict(rid{i}).gene;
    lmIdx(i) = dict(rid{i}).islm;
end
%%
for i = 21:numel(batches)
    batch = batches{i};
    fprintf('%s %d\n',batch,i);
    filter = BasicDBObject();
    filter.append('det_plate',Pattern.compile(batch));
    plates = readColl.distinct('det_plate',filter);
    plates = j2m(plates);
    
%     % replcate data on this plate has been removed from db
%     idx = strcmp('LJP005_MCF10A_3H_X2_B17',plates);
%     plates(idx) = [];
    platesRes = cell(1,numel(plates));
    
    tic
    % compute chdir for replicates on each plate
    for j = 1:numel(plates)
        plate = plates{j};
        query = BasicDBObject();
        query.append('det_plate',plate);
        cursor = readColl.find(query);
        arr = cell(cursor.count,1);
        cids = cell(cursor.count,1);
        for k = 1:cursor.count
            arr{k} = j2m(cursor.next());
            cids{k} = arr{k}.distil_id;
        end
        t = parse_gctx(bigMatPath,'cid',cids);
        % t.cid matches the order of cids
        for k = 1:numel(cids)
            arr{k}.data = t.mat(:,k);
        end
     
        plateRes = getChdir_2(arr,lmIdx);
        platesRes{j} = plateRes';
    end
    toc
    
    % get unique experiments' sig_ids.
    query = BasicDBObject();
    query.append('sig_id',Pattern.compile(batch));
    query.append('brew_prefix',batch);
    sub = BasicDBObject();
    sub.append('$ne','ctl_vehicle');
    query.append('pert_type',sub);
    jsonProjection = sprintf('{_id:0,sig_id:1,distil_id:1,pert_mfc_id:1}');
    cursor = sigColl.find(query,JSON.parse(jsonProjection));
    sigIdStructs = cell(cursor.count,1);
    for sigIdx = 1:cursor.count
        sigIdStructs{sigIdx} = j2m(cursor.next());
    end
    
    % merge chdir replicates and computeSignficance
    chdirArr = mergeReplicates(platesRes,sigIdStructs);
    
    % save chdir arr to db.
    tic
    for j = 1:numel(chdirArr)
        chdirStruct = chdirArr{j};
        writeColl.save(json2dbobj(savejson('',chdirStruct)));
    end
    toc
end
