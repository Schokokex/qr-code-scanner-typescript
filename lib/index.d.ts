export declare class QRScanner{
    public static initiate: (opts: {
        onResult: (result:any)=>any;
        onError?: (err:any)=>any;
        onTimeout?: ()=>any;
        match?: string;
        timeout?: number;
        parent?: HTMLElement;
        lockLayerParent?: HTMLElement;
        className?: string;
        lockLayerClassName?: string;
    }) => any;
}
